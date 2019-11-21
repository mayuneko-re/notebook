# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
# from numba import jit

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

    R = 8.3142 # Gas constant, kPa-m3/(Kg-mol-K)

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


    def get_Vm(self, P, T):
        """Calculate molar volume at P and T using Peng-Robinson EOS

        Args:
          P (float): Pressure condition, BarA
          T (float): Temperature condition, degK

        Returns:
          float: Molar volume Vm, m3/mol
        """
        z = self.get_z_factor(P, T)
        Vm = z * self.R * T / (P*1e5) # m3/mol
        return Vm



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

    R = 8.3142 # Gas constant, kPa-m3/(Kg-mol-K)

    def __init__(self, feed):
        self.feed = np.array(feed) # Feed component instance list
        self.set_BIPs() # Set zero BIPs


    def set_composition(self, Zi):
        """Set or change composition in feed components
        """
        self.Zi = Zi / np.sum(Zi) # Feed molar fraction
        self.Zi[self.Zi<1e-4] = 0
        self.Zi = self.Zi / np.sum(self.Zi)
        self.Nc = self.feed.size # No of component
        
        if self.feed[self.Zi>0].size == 1:
            self.isPure = True
            self.feed_pure = self.feed[self.Zi>0][0]
        else:
            self.isPure = False



    def set_BIPs(self, kik=None):
        """Set Binary Interaction Parameters for PR EOS. Available from literature.
        """
        if kik is None:
            self.kik = np.zeros((self.feed.size, self.feed.size))
        else:
            self.kik = kik


    def PT_flash(self, P, T, initialize_k=True, verbose=False):
        """PT Flash calculation
        """

        # Pressure and Temperature condition
        self.P = P # BarA
        self.T = T # degK

        # Single component -> not mixture
        if self.isPure:
            # print(self.feed_pure)
            fp = self.feed_pure
            # Geuss Ki
            Ki = (fp.Pc/self.P)*np.exp(5.37*(1+fp.omega)*(1-fp.Tc/self.T))
            zM = fp.get_z_factor(self.P, self.T)[0]
            Vm = fp.get_Vm(self.P, self.T)
            zL, zV = np.nan, np.nan
            self.Xi, self.Yi = np.nan, np.nan
            
            if Ki>1:    
                self.V, self.L = np.inf, -np.inf
            else:
                self.V, self.L = -np.inf, np.inf

            self.phase, self.Vact, self.Lact, self.Xiact, self.Yiact, self.zLact, self.zVact = \
                self.set_output_actual(Ki, self.Xi, self.Yi, self.V, self.L, zM, zL, zV)
            self.Ki = np.nan
 
            self.VmVact = self.zVact * self.R * self.T / (self.P*1e5) # m3/mol
            self.VmLact = self.zLact * self.R * self.T / (self.P*1e5) # m3/mol

            if verbose:
                self.print_result()
            return

        # Get components actually exist        
        Ncm = self.feed[self.Zi>0].size

        if Ncm == self.Nc:
            # PT Flash
            Ki, Xi, Yi, V, L, zV, zL, zM = self.PT_flash_core(self.feed, self.Zi, self.kik, initialize_k=initialize_k)
        else:
            mask = self.Zi>0
            feedm = self.feed[mask]
            Zm = self.Zi[mask]
            Ncm = self.feed[mask].size
            kikm = self.kik[np.outer(mask,mask)].reshape(Ncm,Ncm)       
        
            # PT Flash
            Kim, Xim, Yim, V, L, zV, zL, zM = self.PT_flash_core(feedm, Zm, kikm, initialize_k=initialize_k)
        
            # loc = np.where(self.Zi == 0)[0]
            # Xi = np.insert(Xi, loc, 0)
            # Yi = np.insert(Yi, loc, 0)
            # Ki = np.insert(Ki, loc, np.nan)
            loc = np.where(self.Zi != 0)[0]
            Xi, Yi = np.zeros_like(self.Zi), np.zeros_like(self.Zi)
            Ki = np.zeros_like(self.Zi) * np.nan
            Xi[loc] = Xim
            Yi[loc] = Yim
            Ki[loc] = Kim
        
        # print('z factor of mixture if one phase',zM)

        # Set outputs
        self.Ki =  Ki
        self.Xi, self.L, self.zL = Xi, L, zL
        self.Yi, self.V, self.zV = Yi, V, zV        
        self.VmV = zV * self.R * self.T / (self.P*1e5) # m3/mol
        self.VmL = zL * self.R * self.T / (self.P*1e5) # m3/mol
        self.phase, self.Vact, self.Lact, self.Xiact, self.Yiact, self.zLact, self.zVact = \
            self.set_output_actual(Ki, Xi, Yi, V, L, zM, zL, zV)
        self.VmVact = self.zVact * self.R * self.T / (self.P*1e5) # m3/mol
        self.VmLact = self.zLact * self.R * self.T / (self.P*1e5) # m3/mol

        # Print results
        if verbose:
            self.print_result()
            
        return True


    
    def PT_flash_core(self, feed, Zi, kik, initialize_k=True):
        """PT Flash core engine only for actually exsiting components
        
        Args:
          feed:
          Zi:
          kik:

        Returns:
          Ki: Equilibrium constants 
          Xi: Mole fraction in liquid phase
          Yi: Mole fraction in vapor phase
          V: 
          L:
        """

        # Initialize PR EOS parameters for each components at P and T condition
        for c in feed:
            c.PREOS(self.P, self.T)

        # Mixture parameters are calculated by mixing rules.
        Ai = np.array([c.A for c in feed])
        # print(Ai)
        Bi = np.array([c.B for c in feed])
        Aik = np.outer(Ai,Ai)**0.5 * (1-kik)
        # print(Aik)

        # New values of Ki thus calculated are again used to estimate V and 
        # thereafter Xi & Yi. Iteration is repeated till there is no further 
        # change in Ki values.

        if initialize_k is True:
            # Initial guess of Ki is made by Wilson equation.
            Ki = np.array([(c.Pc/self.P)*np.exp(5.37*(1+c.omega)*(1-c.Tc/self.T)) for c in feed])
        else:
            # Otherwise, set Ki explicitly
            Ki = initialize_k

        def f(V, Z, Ki):
            f = np.sum( (Z*(Ki-1)) / (1+V*(Ki-1)) )
            return f

        def df(V, Z, Ki):
            # For Newton-Raphson method
            df = -np.sum( (Z*(Ki-1)**2) / (1+V*(Ki-1))**2 )
            return df

        deltaKi = 10
        tol = 1e-6

        V = 0.5

        while deltaKi>tol:

            # Single phase
            if np.all(Ki<1):
                # self.phase = 'liquid'
                L, V = np.inf, -np.inf
                Xi, Yi = np.nan, np.nan
                break
            elif np.all(Ki>1):
                # self.phase = 'vapor'
                L, V = -np.inf, np.inf
                Xi, Yi = np.nan, np.nan
                break

            # Solve Relative molar volume in vapor phase
            min = 1/(1-np.max(Ki))
            max = 1/(1-np.min(Ki))
            min = min + np.abs(min)*1e-6
            max = max - np.abs(max)*1e-6
            V = optimize.bisect(f, min, max, args=(Zi, Ki)) # Bisection method for Vapor 
            # V = optimize.newton(f, 0.5, df, args=(Zi, Ki)) # Newton-Raphson method for Vapor 
            # V = optimize.newton(f, V, df, args=(Zi, Ki)) # Newton-Raphson method for Vapor 
            L = 1 - V # Liquid

            # Compositions in liquid phase and vapor phase
            Xi = Zi / (1+V*(Ki-1)) # Liquid
            Yi = Xi * Ki # Vapor

            # Partial fugacity coefficient calculation for Liquid phase and Vapor phase
            PhiL, zL = calc_fugacity_coefficient(Aik, Bi, Xi, self.P, 'liquid')
            PhiV, zV = calc_fugacity_coefficient(Aik, Bi, Yi, self.P, 'vapor')

            # New K
            Kinew = PhiL / PhiV
            deltaKi = np.sum(np.abs(Kinew/Ki-1))
            Ki = Kinew

        if V>1:
            _, zM = calc_fugacity_coefficient(Aik, Bi, Zi, self.P, 'vapor')
        elif V<0:
            _, zM = calc_fugacity_coefficient(Aik, Bi, Zi, self.P, 'liquid')
        else:
            zM = np.nan

        return Ki, Xi, Yi, V, L, zV, zL, zM





    def set_output_actual(self, Ki, Xi, Yi, V, L, zM, zL, zV):
        
        # res = [{phase: [Vact, Lact, Xiact, Yiact, zLact, zVact], ...}]
        res = {
            'vapor':    [1, 0, np.nan , self.Zi, np.nan, zM    ],
            'liquid':   [0, 1, self.Zi, np.nan , zM    , np.nan],
            'two-phase':[V, L, Xi     , Yi     , zL    , zV    ]
            }

        if V>=1:
            out = 'vapor', *res.get('vapor')
        elif V<=0:
            out = 'liquid', *res.get('liquid')
        else:
            out = 'two-phase', *res.get('two-phase')

        return out


    def print_result(self):
        """Print results.
        """

        # print('PT Flash calculation converged.')
        print()
        print('                       Feed components : {}'.format([c.name for c in self.feed]))
        print('                    Feed mole fraction : {}'.format(self.Zi))
        print()
        print('PT Flash at {0:.1f} degK and {1:.1f} BarA.'.format(self.T, self.P))
        print()
        print('                                 Phase : {}'.format(self.phase))
        print('Relative mole fraction of liquid phase : {:.4f}'.format(self.Lact))
        print(' Relative mole fraction of vapor phase : {:.4f}'.format(self.Vact))
        print('              z factor of liquid phase : {:.4f}'.format(self.zLact))
        print('               z factor of vapor phase : {:.4f}'.format(self.zVact))
        print('          Molar volume of liquid phase : {:.6f} m3/mol'.format(self.VmLact))
        print('           Molar volume of vapor phase : {:.6f} m3/mol'.format(self.VmVact))
        print('         Mole fraction in liquid phase : {}'.format(self.Xiact))
        print('          Mole fraction in vapor phase : {}'.format(self.Yiact))
        print('                              K values : {}'.format(self.Ki))

        return True



def calc_fugacity_coefficient(Aik, Bi, xi, P, phase):
    # phase: 'vapor' or 'liquid'

    # Mixture parameters are calculated by mixing rules.
    A = np.sum(Aik * xi * xi.reshape(-1, 1))

    # print(xi,A)

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

    return Phi, Zj


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
