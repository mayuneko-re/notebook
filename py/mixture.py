# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


class Components():

    # Gas constant
    R = 8.3142 # kPa-m3/(Kg-mol-K)

    def __init__(self, name, Pc, Tc, omega):
        self.name = name
        self.Pc = Pc # BarA
        self.Tc = Tc # degK
        self.omega = omega # acentric factor

    def PREOS(self, P, T):
        # Calculate Peng-Robinson EOS parameters
        
        self.P = P # BarA
        self.T = T # degK

        m = 0.37464 + 1.54226*self.omega -0.26992*self.omega**2
        self.alpha = (1 + m*(1-np.sqrt(T/self.Tc)))**2
        self.a = 0.45724*(self.R*self.Tc)**2*self.alpha/(self.Pc*1e5)
        self.b = 0.07780*self.R*self.Tc/(self.Pc*1e5)
        self.A = self.a * (self.P*1e5) / (self.R*self.T)**2
        self.B = self.b * (self.P*1e5) / (self.R*self.T)
        
        return True



class Mixture():

    # Gas constant
    R = 8.3142 # kPa-m3/(Kg-mol-K)

    def __init__(self, feed, Zi):
        self.feed = np.array(feed) # Feed composition instance list
        self.Zi = Zi / np.sum(Zi) # Feed molar fraction
        self.set_BIPs() # Set zero BIPs

    def set_BIPs(self, kik=None):
        # Set Binary Interaction Parameters for PR EOS.
        # Available from literature.
        if kik is None:
            self.kik = np.zeros((self.feed.size, self.feed.size))
        else:
            self.kik = kik
            
            
    def PT_flash(self, P, T, verbose=False):
        # PT Flash

        self.P = P # BarA
        self.T = T # degK

        # Get components actually exist
        self.feedm = self.feed[self.Zi>0]
        self.Zm = self.Zi[self.Zi>0]
        self.Ncm = self.feed[self.Zi>0].size
        # print(self.kik)
        self.kik = self.kik[np.outer(self.Zi>0,self.Zi>0)].reshape(self.Ncm,self.Ncm)       
        # print(self.kik)
        
        # Single component system
        if self.Ncm == 1:
            self.phase = 'single-component-system'
            self.Vact, self.V = np.nan, np.nan
            self.Lact, self.L = np.nan, np.nan
            self.Xiact, self.Xi = np.nan, np.nan
            self.Yiact, self.Yi = np.nan, np.nan
            self.Ki = np.nan
            if verbose:
                self.print_result()
            return


        # Initialize PR EOS parameters for each components at P and T condition
        for c in self.feedm:
            c.PREOS(self.P, self.T)

        # Initial guess of Ki is made by Wilson equation.
        self.Ki = np.array([(c.Pc/self.P)*np.exp(5.37*(1+c.omega)*(1-c.Tc/self.T)) 
                            for c in self.feedm])

        # Mixture parameters are calculated by mixing rules.
        self.Ai = np.array([c.A for c in self.feedm])
        self.Bi = np.array([c.B for c in self.feedm])
        self.Aik = np.outer(self.Ai,self.Ai)**0.5 * (1-self.kik)
        
   
        def f(V, Z, Ki):
            f = np.sum( (Z*(Ki-1)) / (1+V*(Ki-1)) )
            return f

        # New values of Ki thus calculated are again used to estimate V and 
        # thereafter Xi & Yi. Iteration is repeated till there is no further 
        # change in Ki values.
        deltaKi = 10
        tol = 1e-6

        while deltaKi>tol:

            # Single component system
            if np.all(self.Ki<1):
                self.phase = 'liquid'
                self.Vact, self.V = 0, np.nan
                self.Lact, self.L = 1, np.nan
                self.Xiact, self.Xi = self.Zi, np.nan
                self.Yiact, self.Xi = np.nan, np.nan
                self.Ki = np.nan
                if verbose:
                    self.print_result()
                return

            elif np.all(self.Ki>1):
                self.phase = 'vapor'
                self.Vact, self.V = 1, np.nan
                self.Lact, self.L = 0, np.nan
                self.Xiact, self.Xi = np.nan, np.nan
                self.Yiact, self.Yi = self.Zi, np.nan
                self.Ki = np.nan
                if verbose:
                    self.print_result()
                return
            
            # Relative molar volume in vapor phase
            # self.V = optimize.fsolve(f, 0.5, args=(self.Zi, self.Ki))
            min = 1/(1-np.max(self.Ki))
            max = 1/(1-np.min(self.Ki))
            min = min + np.abs(min)*1e-6
            max = max - np.abs(max)*1e-6

            self.V = optimize.bisect(f, min, max, args=(self.Zm, self.Ki))
            self.L = 1 - self.V

            # Compositions in liquid phase and vapor phase
            self.Xi = self.Zm / (1+self.V*(self.Ki-1))
            self.Yi = self.Xi * self.Ki

            # Partial fugacity coefficient calculation for Liquid phase and Vapor phase
            PhiL = self.calc_fugacity(self.Xi, P)
            PhiV = self.calc_fugacity(self.Yi, P)

            # New K
            Kinew = PhiL / PhiV
            deltaKi = np.sum(np.abs(Kinew/self.Ki-1))
            self.Ki = Kinew

        loc = np.where(self.Zi == 0)[0]
        self.Xi = np.insert(self.Xi, loc, 0)
        self.Yi = np.insert(self.Yi, loc, 0)
        self.Ki = np.insert(self.Ki, loc, np.nan)

        if self.V >=1:
            self.phase = 'vapor'
            self.Vact = 1
            self.Lact = 0
            self.Xiact = np.zeros_like(self.Zi)
            self.Yiact = self.Zi
        elif self.V<=0:
            self.phase = 'liquid'
            self.Vact = 0
            self.Lact = 1
            self.Xiact = self.Zi
            self.Yiact = np.zeros_like(self.Zi)
        else:
            self.phase = 'two-phase'
            self.Vact = self.V
            self.Lact = self.L
            self.Xiact = self.Xi
            self.Yiact = self.Yi

        if verbose:
            self.print_result()
            



    def calc_fugacity(self, xi, P):

        # Mixture parameters are calculated by mixing rules.
        A = np.sum(self.Aik * xi * xi.reshape(-1, 1))
        B = np.sum(self.Bi * xi)

        Zj = CardanoEOS(A,B)

        lnPhi = self.Bi/B*(Zj-1) - np.log(Zj-B) \
            - A/(2*np.sqrt(2)*B)*(2*np.dot(self.Aik, xi)/A-self.Bi/B) \
            * np.log( (Zj+(1+np.sqrt(2))*B)/(Zj+(1-np.sqrt(2))*B) )
        Phi = np.exp(lnPhi)

        return Phi


    def print_result(self):
        print('PT Flash calculation converged.')
        print()
        print('                       Feed components : {}'.format([c.name for c in self.feed]))
        print('                    Feed mole fraction : {}'.format(self.Zi))
        print()
        print('PT Flash at {0:.1f} degK and {1:.1f} BarA.'.format(self.T, self.P))
        print()
        print('                                 Phase : {}'.format(self.phase))
        print('Relative mole fraction of liquid phase : {:.4f}'.format(self.Lact))
        print(' Relative mole fraction of vapor phase : {:.4f}'.format(self.Vact))
        print('         Mole fraction in liquid phase : {}'.format(self.Xiact))
        print('          Mole fraction in vapor phase : {}'.format(self.Yiact))
        print('                              K values : {}'.format(self.Ki))

def CardanoEOS(A,B):

    C2 = B-1
    C1 = A - 3*B**2 - 2*B
    C0 = B**3 + B**2 - A*B
    q = C0 - 1/3*C1*C2 + 2/27*C2**3
    p = C1 - C2**2/3
    D = (q/2)**2 + (p/3)**3

    return np.cbrt(-q/2+np.sqrt(D)) + np.cbrt(-q/2-np.sqrt(D)) - C2/3

