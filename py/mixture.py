# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# ------------------------#
#      Componnents        #
# ------------------------#
class Components():
    """Set component.

    Attributes:
        name (str): Name of component
        Pc (float): Critical pressure, BarA
        Tc (float): Critical temperature, degK
        omega (float): Acentric factor
        P (float): Pressure condition, BarA
        T (float): Temperature condition, degK
        alpha (float): PR EOS parameter at T
        a (float): PR EOS parameter at T
        b (float): PR EOS parameter
        A (float): PR EOS parameter at P and T
        B (float): PR EOS parameter at P and T
    """

    # Gas constant
    R = 8.3142 # kPa-m3/(Kg-mol-K)

    def __init__(self, name, Pc, Tc, omega):
        self.name = name
        self.Pc = Pc # BarA
        self.Tc = Tc # degK
        self.omega = omega # acentric factor


    def PREOS(self, P, T):
        """Calculate Peng-Robinson EOS parameters at P and T.
        
        Args:
          P (float): Pressure condition, BarA
          T (float): Temperature condition, degK
        """
        
        self.P = P # BarA
        self.T = T # degK

        m = 0.37464 + 1.54226*self.omega - 0.26992*self.omega**2
        self.alpha = (1 + m*(1-np.sqrt(T/self.Tc)))**2
        self.a = 0.45724*(self.R*self.Tc)**2*self.alpha/(self.Pc*1e5)
        self.b = 0.07780*self.R*self.Tc/(self.Pc*1e5)
        self.A = self.a * (self.P*1e5) / (self.R*self.T)**2
        self.B = self.b * (self.P*1e5) / (self.R*self.T)
        
        return True

    def get_z_factor(self, P, T):
        """Calculate z factor at P and T using Peng-Robinson EOS .
        
        Args:
          P (float): Pressure condition, BarA
          T (float): Temperature condition, degK

        Returns:
          float: z factor (one or three real values)
        """

        self.P = P # BarA
        self.T = T # degK

        self.PREOS(P,T)

        z = find_z_factor(self.A,self.B)

        """
        if phase == 'vapor':
            Zj = np.max(Zj)
        elif phase == 'liquid':
            Zj = np.min(Zj)
        """

        return z





# ------------------------#
#        Mixture          #
# ------------------------#
class Mixture():
    """Set mixture.

    Attributes:
        feed (ndarray): ndarray of feed component list (Instances of Components Class)
        Zi (ndarray): Feed molar fraction for each feed component
        kik (ndarray): Binary interaction parameters matrix for Peng-Robinson EOS
        P (float): Pressure condition for PT Flash, BarA
        T (float): Temperature condition for PT Flash, degK
    """

    # Gas constant
    R = 8.3142 # kPa-m3/(Kg-mol-K)

    def __init__(self, feed):
        self.feed = np.array(feed) # Feed component instance list
        self.set_BIPs() # Set zero BIPs


    def set_composition(self, Zi):
        """Set or change composition in feed components
        """
        self.Zi = Zi / np.sum(Zi) # Feed molar fraction


    def set_BIPs(self, kik=None):
        """Set Binary Interaction Parameters for PR EOS. Available from literature.
        """
        if kik is None:
            self.kik = np.zeros((self.feed.size, self.feed.size))
        else:
            self.kik = kik
            
            
    def PT_flash(self, P, T, verbose=False):
        # PT Flash

        # Pressure and Temperature condition
        self.P = P # BarA
        self.T = T # degK

        # Get components actually exist
        self.feedm = self.feed[self.Zi>0]
        self.Zm = self.Zi[self.Zi>0]
        self.Ncm = self.feed[self.Zi>0].size
        self.kikm = self.kik[np.outer(self.Zi>0,self.Zi>0)].reshape(self.Ncm,self.Ncm)       
        
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
        self.Aik = np.outer(self.Ai,self.Ai)**0.5 * (1-self.kikm)
        
   
        def f(V, Z, Ki):
            f = np.sum( (Z*(Ki-1)) / (1+V*(Ki-1)) )
            return f

        # New values of Ki thus calculated are again used to estimate V and 
        # thereafter Xi & Yi. Iteration is repeated till there is no further 
        # change in Ki values.
        deltaKi = 10
        tol = 1e-6

        while deltaKi>tol:

            # Single phase
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
            
            # Solve Relative molar volume in vapor phase
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
            PhiL = calc_fugacity_coefficient(self.Aik, self.Bi, self.Xi, P, 'vapor')
            PhiV = calc_fugacity_coefficient(self.Aik, self.Bi,self.Yi, P, 'liquid')

            # New K
            Kinew = PhiL / PhiV
            deltaKi = np.sum(np.abs(Kinew/self.Ki-1))
            self.Ki = Kinew


        loc = np.where(self.Zi == 0)[0]
        self.Xi = np.insert(self.Xi, loc, 0)
        self.Yi = np.insert(self.Yi, loc, 0)
        self.Ki = np.insert(self.Ki, loc, np.nan)

        # Set outputs
        # {phase: [Vact, Lact, Xiact, Yiact]}
        res = {
            'vapor':    [1     , 0     , np.zeros_like(self.Zi), self.Zi               ],
            'liquid':   [0     , 1     , self.Zi               , np.zeros_like(self.Zi)],
            'two-phase':[self.V, self.L, self.Xi               , self.Yi               ]
            }

        if self.V>=1:
            out = 'vapor', *res.get('vapor')
        elif self.V<=0:
            out = 'liquid', *res.get('liquid')
        else:
            out = 'two-phase', *res.get('two-phase')

        self.phase, self.Vact, self.Lact, self.Xiact, self.Yiact = out

        # Print results
        if verbose:
            self.print_result()
            
        return True












    def print_result(self):
        """Print results.
        """

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

        return True





def calc_fugacity_coefficient(Aik, Bi, xi, P, phase):
    # phase: 'vapor' or 'liquid'

    # Mixture parameters are calculated by mixing rules.
    A = np.sum(Aik * xi * xi.reshape(-1, 1))
    B = np.sum(Bi * xi)

    Zj = find_z_factor(A,B)
    if phase == 'vapor':
        Zj = np.max(Zj)
    elif phase == 'liquid':
        Zj = np.min(Zj)


    lnPhi = Bi/B*(Zj-1) - np.log(Zj-B) \
        - A/(2*np.sqrt(2)*B)*(2*np.dot(Aik, xi)/A-Bi/B) \
        * np.log( (Zj+(1+np.sqrt(2))*B)/(Zj+(1-np.sqrt(2))*B) )
    Phi = np.exp(lnPhi)

    return Phi






def find_z_factor(A,B):
    """Solves a Cubic equation for z factor.
    z^3 + (B-1) z^2 + (A-3*B**2-2*B) z + (B**3+B**2-A*B) = 0 

    Args:
      A (float): See above.
      B (float): See above.

    Returns:
      float: z factor (one or three real values)
    """

    C2 = B - 1
    C1 = A - 3*B**2 - 2*B
    C0 = B**3 + B**2 - A*B
    
    # z = Cardano(C2,C1,C0)
    z = np.roots([1,C2,C1,C0])
    z = z[np.isreal(z)].real

    return z









def Cardano(C2,C1,C0):
    """Solves a Cubic equation by Cardano method
    x^3 + C2 x^2 + C1 x + C0 = 0

    Args:
      C2 (float): See above.
      C1 (float): See above.
      C0 (float): See above.

    Returns:
      float: Root for a Cubic equation
    """

    q = C0 - 1/3*C1*C2 + 2/27*C2**3
    p = C1 - C2**2/3
    D = (q/2)**2 + (p/3)**3
    x = np.cbrt(-q/2+np.sqrt(D)) + np.cbrt(-q/2-np.sqrt(D)) - C2/3

    return x
