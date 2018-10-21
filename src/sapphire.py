#!/usr/local/bin/python

#python Version 2.7.2
import numpy as np
import glob
import pickle as pkl
import scipy.integrate as itg
import copy

#Class for sapphire
class Sapphire:
    def __init__(self):
        #***** Private objects *****
        self.__ph = Physics()

        #***** Public variables *****
        #Indexes of refraction (measurement)
        self.oN   = 3.05
        self.eN   = 3.38
        self.avgN = (self.oN + self.eN)/2.
        self.oN_err   = 0.03
        self.eN_err   = 0.03
        self.avgN_err = self.__ph.invVar([self.oN_err, self.eN_err]) 
        #Loss tangent (measurement)
        self.oLT   = 0.1e-4
        self.eLT   = 1.1e-4
        self.avgLT = (self.oLT + self.eLT)/2.
        self.oLT_err   = 1.3e-4
        self.eLT_err   = 1.3e-4
        self.avgLT_err = self.__ph.invVar([self.oLT_err, self.eLT_err])
        #Thickness (measurement) [m]
        self.Th = 3.75e-3
        self.Th_err = 0.01e-3

        return
    
    #Sapphire loss tangent vs temperature
    def ltVsTemp(self, tempArr):
        #Second-order polynomial parameters determined from a best fit of V. V. Parshin
        p = [2.86925077e-09, -8.02203409e-08, -8.51669651e-06]
        #First-order polynomial parameters determined from a best fit of V. V. Parshin
        #p = [8.38487794e-07, -7.10412123e-05]
        poly = np.poly1d(p)
        LT = poly(tempArr)
        if LT < 1.0e-15:
            LT = 0.0
        return LT

    #Sapphire loss tangent vs frequency
    def ltVsFreq(self, freq):
        #First-order polynomial fit
        pO = [0.0240*1.e-4, -1.999*1.e-4]
        pE = [0.0237*1.e-4, -1.063*1.e-4]
        LTO = pO[0]*freq + pO[1]
        LTE = pE[0]*freq + pE[1]
        if LTO < 1.e-15:
            LTO = 0.
        if LTE < 1.e-15:
            LTE = 0.
        return LTO, LTE
