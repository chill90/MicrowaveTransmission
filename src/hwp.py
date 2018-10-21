#!/usr/local/bin/python

#python Version 2.7.2
import numpy as np
import glob
import pickle as pkl
import scipy.integrate as itg
import copy

#Class for the HWP
class HWP:
    def __init__(self):
        #***** Private objects *****
        self.__sp = Sapphire()
        self.__ph = Physics()
        
        #***** Public variables *****
        #Number of plates (PB2)
        self.numPlates = 3
        #Default polarization fraction
        self.defPfrac = 1.
        #Default input polarization angle [deg]
        self.defPolAng = 0.
        #Default HWP angle [deg]
        self.defHwpAng = 0.
        #Ideal frequency for the HWP (PB2)
        self.freqIdeal = 120.e9 #[Hz]
        #Ideal plate thickness
        self.thickIdeal = self.hwpThick(self.freqIdeal)
        #Temperature -- 273 K for PB2a, 100 K for PB2b
        self.warmTemp = 273.
        self.coldTemp = 100.
        #Ideal plate angles [deg]
        self.plateAngles_1 = np.array([0.])
        self.plateAngles_3 = np.array([0., 54., 0.])
        self.plateAngles_5 = np.array([0., 29., 94.5, 29., 2.])
        #Ideal plate interface angles
        self.interfaceAngles_1 = np.array([0.])
        self.interfaceAngles_3 = np.array([self.plateAngles_3[i+1]-self.plateAngles_3[i] for i in range(len(self.plateAngles_3)-1)])
        self.interfaceAngles_5 = np.array([self.plateAngles_5[i+1]-self.plateAngles_5[i] for i in range(len(self.plateAngles_5)-1)])
        #Mueller polarizers
        self.__GX = 0.5*np.matrix([[1, 1,0,0], [1, 1,0,0], [0,0,0,0], [0,0,0,0]])
        self.__GY = 0.5*np.matrix([[1,-1,0,0], [-1,1,0,0], [0,0,0,0], [0,0,0,0]])
        
        return
    
    #Mueller matrix formalism using code from Tomo's thesis
    #***** Private functions *****
    #Rotation matrix
    def __R(self, theta):
        theta = self.__ph.degToRad(theta)
        return np.matrix([[1,0,0,0],[0,np.cos(2*theta),-np.sin(2*theta),0],[0,np.sin(2*theta),np.cos(2*theta),0],[0,0,0,1]])

    #Retardance matrix
    def __Gamma(self, delta):
        delta = self.__ph.degToRad(delta)
        return np.matrix([[1,0,0,0],[0,1,0,0],[0,0,np.cos(delta),-np.sin(delta)],[0,0,np.sin(delta),np.cos(delta)]])

    #***** Public Functions *****
    #Half-wave plate thickness [m]
    def hwpThick(self, freq=None, oN=None, eN=None):
        if freq == None:
            freq = self.freqIdeal
        if oN == None:
            oN = self.__sp.oN
        if eN == None:
            eN = self.__sp.eN

        return np.pi/((2*np.pi)*(eN - oN)/self.__ph.lamb(freq))
    
    #Analytic function for calculating modulation phase
    def phaseAnalytic(self, intensity, modEff, alpha, rho):
        return -0.25*self.__ph.radToDeg(np.arccos((2.*intensity - 1.)/(modEff))) + rho - 0.5*alpha
    
    #HWP transformation matrix for a perfect HWP (no internal reflections)
    def T(self, nu, rho, angles=None, d=None, oN=None, eN=None):
        if angles == None:
            angles = self.plateAngles_3
        if d == None:
            d = self.thickIdeal
        if oN == None:
            oN = self.__sp.oN
        if eN == None:
            eN = self.__sp.eN        

        Tret = np.matrix([[1,0,0,0], [0,1,0,0],[0,0,1,0],[0,0,0,1]])
        for i in range(len(angles)):
            Tret = self.__R(-angles[i]-rho)*self.__Gamma(self.__ph.birefringentRot(nu, d[i], oN, eN))*self.__R(angles[i]+rho)*Tret
        return Tret

    #Output Stokes vector for a perfect HWP (no internal reflections)
    def Sout(self, alpha, nu, rho, angles=None, d=None, P=None):
        if angles == None:
            angles = self.plateAngles_3
        if d == None:
            d = self.thickIdeal
        if P == None:
            P = self.defPfrac

        return self.T(nu, rho, angles, d)*self.__ph.Stokes(P,alpha)

    #Modulation efficiency for a perfect HWP (no internal reflections)
    def modEff(self, alpha, nu, angles=None, d=None, P=None):
        if angles == None:
            angles = self.plateAngles_3
        if d == None:
            d = self.thickIdeal
        if P == None:
            P = self.defPfrac

        rho = np.linspace(0., 90., 180)
        Intensities = []
        for i in rho:
            Intensities.append((self.__GX*self.Sout(alpha, nu, i, angles, d, P)).item(0))
        Intensities = np.array(Intensities)
        modEff = (np.amax(Intensities) - np.amin(Intensities))/(np.amax(Intensities) + np.amin(Intensities))
        return modEff

    #Modulation phase for a perfect HWP (no internal reflections) [deg]
    def phase(self, alpha, nu, rho, angles=None, d=None, P=None):
        if plateAngles == None:
            plateAngles = self.plateAngles_3
        if d == None:
            d = self.thickIdeal
        if P == None:
            P = self.defPfrac
        
        intensityRef = (self.__GX*self.Sout(alpha, freq, rho, angles, d, P)).item(0)
        intensity    = (self.__GX*self.Sout(alpha, nu,   rho, angles, d, P)).item(0)
        phaseRef = -0.25*self.__ph.radToDeg(np.arccos((2.*intensityRef - 1.)/(MEperfect(alpha, freq, angles, d, P)))) + rho - 0.5*alpha
        phase    = -0.25*self.__ph.radToDeg(np.arccos((2.*intensity    - 1.)/(MEperfect(alpha, nu,   angles, d, P)))) + rho - 0.5*alpha 
        return phase - phaseRef
