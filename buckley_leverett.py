# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

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
        vD (ndarray): Wave velocity
        vD_sf (float): Shock Front velocity
        Sw_sf (float): Shock Front saturation
        fw_wSF (ndarray): Fractional flow of water with Shock Front
        tD (ndarray): Dimensionless time for plot
        tD_BT (float): Dimensionless time at breakthrough
        Sw_out (ndarray): Water saturation at outlet
        fw_out (ndarray): Fractional flow of water at outlet
        N (ndarray): Oil production volume, PV
        N_BT (float): Oil production volume at breakthrough, PV
        RF (ndarray): Oil recovery factor, HCPV
        RF_BT (float): Oil recovery factor, HCPV
    """
    
    def __init__(self, muw, muo, Swc, Sor, krw0, kro0, nw, no, label=None):

        self.muw = muw # cP, Water viscosity
        self.muo = muo # cP, Oil viscosity
        self.Swc = Swc # Vol/Vol, Connate water saturation
        self.Sor = Sor # Vol/Vol, Residual oil saturation
        self.krw0 = krw0 # Water endpoint relative permeability at residual oil saturation
        self.kro0 = kro0 # Oil endpoint relative permeability at connate water saturation
        self.nw = nw # Exponent for water relative permeability
        self.no = no # Exponent for oil relative permeability
        self.label = label # Label for plot

    def calc(self):

        # Water saturation and normalized water saturation
        self.Sw = np.linspace(start=self.Swc, stop=1-self.Sor, num = 300)
        self.Swn = (self.Sw-self.Swc) / (1-self.Swc-self.Sor)

        self.krw = self.__relperm(Sn=self.Swn, kr0=self.krw0, n=self.nw) # Water
        self.kro = self.__relperm(Sn=1-self.Swn, kr0=self.kro0, n=self.no) # Oil

        # Fractional flow of water
        # fw = 1 / (1 + (kro/oil_viscosity) / (krw/water_viscosity))
        self.fw = (self.krw/self.muw) / (self.krw/self.muw + self.kro/self.muo) # To avoid division by zero

        # Wave velocity = Derivative of fractional flow of water
        self.vD = np.gradient(self.fw, self.Sw)

        # Shock Front
        slope_fw = (self.fw - 0) / (self.Sw - self.Swc)  # Warning: slope_fw[0] = nan
        self.vD_sf = np.nanmax(slope_fw) # Shock Front velocity
        self.Sw_sf = self.Sw[np.nanargmax(slope_fw)] # Shock Front water saturation

        # Dimensionless velocity with Shock Front
        self.vD[self.Sw<self.Sw_sf] = self.vD_sf

        # Fractional flow of water with Shock Front
        self.fw_wSF = np.where(self.Sw<=self.Sw_sf, self.vD_sf*(self.Sw-self.Swc), self.fw)

        self.tD = np.linspace(start=0, stop=2, num = 100)
        
        # Dimensionless breakthrough time
        self.tD_BT = 1/self.vD_sf

        # Water saturation at outlet
        f = interpolate.interp1d(self.vD, self.Sw, kind='linear', bounds_error=False)
        self.Sw_out = np.where(self.tD<=self.tD_BT, self.Swc, f(1/self.tD))

        # Fractional flow of water at outlet
        f = interpolate.interp1d(self.vD, self.fw, kind='linear', bounds_error=False)
        self.fw_out = np.where(self.tD<=self.tD_BT, 0, f(1/self.tD))

        # Average water saturation in porous media
        self.Sw_ave = np.where(self.tD<=self.tD_BT, self.Swc+self.tD, self.Sw_out-(self.fw_out-1)/(1/self.tD))

        self.N = np.where(self.tD<=self.tD_BT, self.tD, self.Sw_ave-self.Swc)
        self.N_BT = self.tD_BT

        self.RF = self.N / (1 - self.Swc)
        self.RF_BT = self.N_BT / (1 - self.Swc) 
    
    # Relative permeability model
    def __relperm(self, Sn, kr0, n):
        kr = kr0 * Sn ** n
        return kr

    # Get saturation profile at specific dimensionless time
    def profileSw(self, tD):
        xD = self.vD * tD
        xD = [1, *xD]
        Sw = [self.Swc, *self.Sw]
        return xD, Sw
