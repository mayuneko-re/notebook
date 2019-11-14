# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate

np.seterr(divide='ignore', invalid='ignore')

class BL:
    """Store and calculate data for oil-water displacement problem by Buckley-Leverett Solution
    
    Attributes:
        muw (float): Water viscosity, cP
        muo (float): Oil viscosity, cP
        Swc (float): Connate water saturation, Vol/Vol
        Sor (float): Residual oil saturation, Vol/Vol
        krw0 (float): Water endpoint relative permeability at residual oil saturation
        kro0 (float): Oil endpoint relative permeability at connate water saturation
        nw (float): Exponent for water relative permeability
        no (float): Exponent for oil relative permeability
        label (str): Case label for plot
        Sw (ndarray): Water saturation
        Swn (ndarray): Normalized water saturation
        krw (ndarray): Water relative permeability
        kro (ndarray): Oil relative permeability
        fw (ndarray): Fractional flow of water
        fw_wSF (ndarray): Fractional flow of water with Shock Front
        vD (ndarray): Wave velocity
        vD_sf (float): Shock Front velocity
        Sw_sf (float): Shock Front saturation
        tD_BT (float): Dimensionless time at breakthrough
        N_BT (float): Oil production volume at breakthrough, PV
        RF_BT (float): Oil recovery factor, HCPV
    """
    
    def __init__(self, muw, muo, Swc, Sor, krw0, kro0, nw, no, label=None):
        """Initialize class
        """
        
        self.muw = muw # cP, Water viscosity
        self.muo = muo # cP, Oil viscosity
        self.Swc = Swc # Vol/Vol, Connate water saturation
        self.Sor = Sor # Vol/Vol, Residual oil saturation
        self.krw0 = krw0 # Water endpoint relative permeability at residual oil saturation
        self.kro0 = kro0 # Oil endpoint relative permeability at connate water saturation
        self.nw = nw # Exponent for water relative permeability
        self.no = no # Exponent for oil relative permeability
        self.label = label # Label for plot
        self.with_gravity = False # Flag for gravity effect

    def calc(self):
        """Calculate values by Buckley-Leverett Solution
        """
        
        # Water saturation and normalized water saturation
        self.Sw = np.linspace(start=self.Swc, stop=1-self.Sor, num = 300)
        self.Swn = (self.Sw-self.Swc) / (1-self.Swc-self.Sor)

        # Relative permeability
        self.krw = self.__relperm(Sn=self.Swn, kr0=self.krw0, n=self.nw) # Water
        self.kro = self.__relperm(Sn=1-self.Swn, kr0=self.kro0, n=self.no) # Oil

        # Fractional flow of water
        self.fw = self.calc_fw()

        # Wave velocity = Derivative of fractional flow of water
        self.vD = np.gradient(self.fw, self.Sw)

        # Determine Shock Front velocity and saturation
        slope_fw = (self.fw - 0) / (self.Sw - self.Swc)  # Warning: slope_fw[0] = nan
        self.vD_sf = np.nanmax(slope_fw) # Shock Front velocity
        self.Sw_sf = self.Sw[np.nanargmax(slope_fw)] # Shock Front water saturation

        # Dimensionless velocity with Shock Front
        self.vD[self.Sw<self.Sw_sf] = self.vD_sf
        self.vD[self.fw>1] = 0

        # Fractional flow of water with Shock Front
        self.fw_wSF = np.where(self.Sw<=self.Sw_sf, self.vD_sf*(self.Sw-self.Swc), self.fw)
        self.fw_wSF[self.fw>1] = 1
        
        # At breakthrough
        self.tD_BT = 1/self.vD_sf # Dimensionless time
        self.N_BT = self.tD_BT # Oil production
        self.RF_BT = self.N_BT / (1 - self.Swc) # Oil recovery factor

    def __relperm(self, Sn, kr0, n):
        """Relative permeability model (Corey shape)
        """
        kr = kr0 * Sn ** n
        return kr

    def calc_fw(self):
        """Fractional flow of water
        """
        # fw = 1 / (1 + (kro/oil_viscosity) / (krw/water_viscosity))
        fw = (self.krw/self.muw) / (self.krw/self.muw + self.kro/self.muo) # To avoid division by zero
        return fw

    def get_Sw_outlet(self, tD):
        """Water saturation at outlet
        """
        # f = interpolate.interp1d(self.vD, self.Sw, kind='linear', bounds_error=False)
        f = interpolate.interp1d(self.vD[::-1], self.Sw[::-1], kind='linear', bounds_error=False, assume_sorted=True)
        # f = interpolate.interp1d(bl.vD[::-1], bl.Sw[::-1], kind='linear', bounds_error=False, assume_sorted=True)
        Sw_out = np.where(tD<=self.tD_BT, self.Swc, f(1/tD))
        return Sw_out

    def get_Fw_outlet(self, tD):
        """Fractional flow of water at outlet
        """
        # f = interpolate.interp1d(self.vD, self.fw, kind='linear', bounds_error=False)
        f = interpolate.interp1d(self.vD, self.fw_wSF, kind='linear', bounds_error=False)
        fw_out = np.where(tD<=self.tD_BT, 0, f(1/tD))
        return fw_out

    def get_Sw_ave(self, tD):
        """Average water saturation in porous media
        """
        Sw_out = self.get_Sw_outlet(tD)
        fw_out = self.get_Fw_outlet(tD)
        Sw_ave = np.where(tD<=self.tD_BT, self.Swc+tD, Sw_out-(fw_out-1)/(1/tD))
        return Sw_ave

    def get_oil_production(self, tD):
        """Oil production volume
        """
        Sw_ave = self.get_Sw_ave(tD)
        N = np.where(tD<=self.tD_BT, tD, Sw_ave-self.Swc)
        return N

    def get_oil_RF(self, tD):
        """Oil recovery factor
        """
        N = self.get_oil_production(tD)
        RF = N / (1 - self.Swc)
        return RF

    def get_Sw_profile(self, tD):
        """Get saturation profile at specific dimensionless time
        """
        xD = self.vD * tD
        xD = [1, *xD]
        Sw = [self.Swc, *self.Sw]
        return xD, Sw
