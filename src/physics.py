#!/usr/local/bin/python

#python Version 2.7.2
import numpy as np
import glob
import pickle as pkl
import scipy.integrate as itg
import copy

#Class for handling general physics
class Physics:
    def __init__(self):
        #***** Private variables
        self.__GHz = 1.e9  #GHz from Hz
        
        # *****Public variables *****
        #Physical Constants
        #Everything is in MKS units
        #Planck constant [J/s]
        self.h = 6.6261e-34
        #Boltzmann constant [J/K]
        self.kB = 1.3806e-23
        #Speed of light [m/s]
        self.c = 299792458.0
        #Pi
        self.PI = 3.14159265
        #Permability of free space [H/m]
        self.mu0 = 1.256637e-6
        #Permittivity of free space [F/m]
        self.ep0 = 8.854188e-12
        #Impedance of free space [Ohm]
        self.Z0 = np.sqrt(self.mu0/self.ep0)
        #CMB Temperature [K]
        self.Tcmb = 2.725
        #Index in vacuum
        self.vacIndex = 1.
        #1/4 phase
        self.quarterPhase = (np.pi/2.)/(2*np.pi)

        #CO Emission lines in [Hz]
        self.coJ10 = 115*self.__GHz
        self.coJ21 = 230*self.__GHz
        self.coJ32 = 345*self.__GHz
        self.coJ43 = 460*self.__GHz

    #***** Public Methods *****
    #Convert from from frequency [Hz] to wavelength [m]
    def lamb(self, freq, index=None):
        if index == None:
            index = self.vacIndex
        
        return self.c/(freq*index)

    #Convert from phase [rad] to physical thickness [m]
    def phaseToThick(self, freq, index=None, phase=None):
        if index == None:
            index = self.vacIndex
        if phase == None:
            phase = self.quarterPhase
        
        return self.lamb(freq, index)*(phase)
    
    #Convert physical thickness [m] to phase [rad]
    def thickToPhase(self, freq, thick=None, index=None):
        if index == None:
            index = self.vacIndex
        if thick == None:
            thick = self.lamb(freq, index)/4.
        
        return (2*np.pi)*(thick/self.lamb(freq, index))

    #Angle rotation by a birefringent material [deg]
    def birefringentRot(self, freq, thick, oN, eN):
        return 360.*(eN - oN)*thick/self.lamb(freq)

    #Stokes vector
    def Stokes(self, polFrac, polAngle):
        polAngle = self.degToRad(polAngle)
        return np.matrix([[1.],[polFrac*np.cos(2.*polAngle)],[polFrac*np.sin(2.*polAngle)],[0]])
        
    #Convert from central frequncy and fractional bandwidth to band edges
    def bandEdges(self, freqCent, fracBw):
        freq1 = freqCent - (freqCent*fracBw)/2
        freq2 = freqCent + (freqCent*fracBw)/2
        return freq1, freq2

    #Spillover efficiency
    def spillEff(self, D, F, waistFact, freq): 
        return 1. - np.exp((-np.power(np.pi,2)/2.)*np.power((D/(waistFact*F*(self.c/freq))),2))

    #Aperture illumination efficiency
    def apertIllum(self, D, F, waistFact, freq):
        lamb = (self.c/freq)
        w0 = D/waistFact
        thetaS = lamb/(np.pi*w0)
        thetaA = np.arctan(1./(2.*F))
        V = lambda theta: np.exp(-np.power(theta, 2.)/np.power(thetaS, 2.))
        effNum = np.power(itg.quad(lambda x: V(x)*np.tan(x/2.), 0., thetaA)[0],2)
        effDenom = itg.quad(lambda x: (np.power(V(x),2.)*np.sin(x)), 0., thetaA)[0]
        effFact = 2.*np.power(np.tan(thetaA/2.),-2.)
        return (effNum/effDenom)*effFact

    #Scattering efficiency off of a rough conductor
    def ruzeEff(self, sigma, f): 
        return np.exp(-np.power(4*np.pi*sigma/(self.c/f), 2.))

    #Ohmic efficiency for reflection off a mirror with finite conductivity
    def ohmicEff(self, sigma, f): 
        return 1. - 4.*np.sqrt(np.pi*f*self.mu0/sigma)/self.Z0

    #Antenna temperature [K] given an intensity [W/m^2] and frequency [Hz]
    def antennaTemp(self, intensity, freq):
        return intensity*(self.c**2)/(2*self.kB*(freq**2))

    #Intensity [W/m^2] given an antenna temperatue [K] and a frequency [Hz]
    def intensityFromAntennaTemp(self, antennaTemp, freq, fbw):
        return 2*(antennaTemp*self.kB*(freq**2)/(self.c**2))*freq*fbw

    #Antenna temperature factor
    def antennaTempFact(self, freq, temp=2.725):
        x = (self.h*freq)/(temp*self.kB)
        thermoFact = np.power((np.exp(x)-1.0),2.0)/(np.power(x,2.0)*np.exp(x))
        return thermoFact

    #Antenna temperature from BB temperature
    def BBToAntennaTemp(self, bbTemp, freq):
        thermoFact = self.antennaTempFact(freq, bbTemp)
        return bbTemp/thermoFact

    #Effective brightness temperature [K]
    def effTemp(self, pow, eff, bc, fbw):
        return pow/(self.kB*bc*fbw*eff)

    #Convert from degrees to radians
    def degToRad(self, deg):
        return (float(deg)/360.)*2*self.PI
    
    #Convert radian to degree
    def radToDeg(self, rad):
        return (rad/(2*self.PI))*360.

    #Inverse variance
    def invVar(self, errArr):
        return 1./np.sqrt(np.sum(1./np.power(np.array(errArr), 2.)))

    #Dielectric loss coefficient with thickness [m] and freq [GHz]
    def dielectricLoss(self, lossTan, thickness, index, freq, atmScatter=0):
        return 1.0 - np.exp((-2*self.PI*index*lossTan*thickness)/self.lamb(freq))

    #Integrate dielectric loss accross a frequency band to obtain emissivity
    def dielectricBandAvgLoss(self, lossTan, thickness, index, freq, fbw):
        freqLo, freqHi = self.bandEdges(freq, fbw) 
        return quad(lambda f: dielectricLoss(lossTan, thickness, index, f), freqLo, freqHi)[0]/(freqHi - freqLo)

    #Rayleigh-Jeans Temperature [K]
    def rjTemp(self, pow, deltaF):
        return pow/(kB*deltaF)

    #Photon Mode Occupation Number
    def nOcc(self, freq, temp):
        return 1./(np.exp((self.h*freq)/(self.kB*temp)) - 1.)

    #Throughput for a diffraction-limited detector [m^2]
    def AOmega(self, freq, nModes=None):
        if nModes == None:
            nModes = 1.

        return nModes*(float(self.c)/freq)**2

    #Blackbody spectral radiance [W/(m^2-Hz)]
    def bbSpecRad(self, emissivity, freq, temp=None):
        if temp == None:
            temp=self.Tcmb

        return emissivity*(2*self.h*(freq**3)/float(self.c**2))*self.nOcc(freq, temp)
    
    #Power spectrum of a blackbody on a diffraction-limited detector [W/Hz]
    def bbPowSpec(self, emissivity, freq, temp=None, nModes=None):
        if temp == None:
            temp=self.Tcmb
        if nModes == None:
            nModes == 1.

        return 0.5*self.AOmega(freq, nModes)*self.bbSpecRad(emissivity, freq, temp)

    #Blackbody power  on a diffraction-limited detector [J]
    def bbPower(self, emissivity, freq, fbw, temp=None, nModes=None):
        if temp == None:
            temp = self.Tcmb
        if nModes == None:
            nModes = 1.

        freq1, freq2 = self.bandEdges(freq, fbw)
        if callable(emissivity):
            return itg.quad(lambda x: self.bbPowSpec(emissivity(x), x, temp, nModes), freq1, freq2)[0]
        else:
            return itg.quad(lambda x: self.bbPowSpec(emissivity   , x, temp, nModes), freq1, freq2)[0]

    #Blackbody power equivalent CMB temperature spectrum on a diffraction-limited detector [K/s]
    def bbPowTempSpec(self, emissivity, freq, temp=None, nModes=None):
        if temp == None:
            temp = self.Tcmb
        if nModes == None:
            nModes = 1.

        return self.bbPowSpec(emissivity, freq, temp)/self.aniPowSpec(emissivity, freq, self.Tcmb)
    
    #Blackbody power equivalent temperature on a diffraction-limited detector [K/s]
    def bbPowTemp(self, emissivity, freq, fbw, temp=None, nModes=None):
        if temp == None:
            temp = self.Tcmb
        if nModes == None:
            nModes = 1.

        freq1, freq2 = bandEdges(freq, fbw)
        if callable(emissivity):
            return itg.quad(self.bbPowTempSpec(emissivity(x), x, temp, nModes), freq1, freq2)[0]
        else:
            return itg.quad(self.bbPowTempSpec(emissivity,    x, temp, nModes), freq1, freq2)[0]

    #Temperature anisotropy power spectrum dP/dT on a diffraction-limited detector [W/K]
    def aniPowSpec(self, emissivity, freq, temp=None):
        if temp == None:
            temp = self.Tcmb

        return ((self.h**2)/self.kB)*emissivity*(self.nOcc(freq, temp)**2)*((freq**2)/(temp**2))*np.exp((self.h*freq)/(self.kB*temp))

    #Temperature anisotropy power dP/dT on a diffraction-limited detector [J/K]
    def aniPower(self, emissivity, freq, fbw, temp=None):
        if temp == None:
            temp = self.Tcmb

        freq1, freq2 = self.bandEdges(freq, fbw)
        if callable(emissivity):
            return itg.quad(lambda x: self.aniPowSpec(emissivity(x), x, temp), freq1, freq2)[0]
        else:
            return itg.quad(lambda x: self.aniPowSpec(emissivity,    x, temp), freq1, freq2)[0]
