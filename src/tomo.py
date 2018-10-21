#!/usr/local/bin/python

#python Version 2.7.2
import numpy as np
import glob
import pickle as pkl
import scipy.integrate as itg
import copy

#Class for an ideal HWP (normal incidence with internal reflections)
#Code taken from Tomo's thesis
class Tomo:
    def __init__(self):
        #Private objects
        self.__hwp = HWP()
        self.__ph  = Physics()
        
        return

    #****** Private functions *****
    def __C(self, nu):
        return 1./(self.__ph.mu0*nu)

    def __Clamb(self, nu, n):
        return self.__C(nu)/self.__ph.lamb(nu, n)

    def __prop(self, nu, n, th):
        return np.exp(1j*self.__ph.thickToPhase(nu, th, n))
        
    def __mprop(self, nu, n, th):
        return np.exp(-1j*self.__ph.thickToPhase(nu, th, n))
    
    def __M1(self, nu, no=None, ne=None):
        if no == None:
            no = self.__hwp.oN
        if ne == None:
            ne = self.__hwp.eN

        Co = __Clamb(nu,no)
        Ce = __Clamb(nu,ne)
        M1Ret = np.matrix([[1,0,1,0],[0,1,0,1],[0,-Ce,0,Ce],[Co,0,-Co,0]])
        return M1Ret
    
    def __M2inv(self, nu, theta, no=None, ne=None, d=None):
        if no == None:
            no = self.__hwp.oN
        if ne == None:
            ne = self.__hwp.eN
        if d == None:
            d = self.__hwp.thickIdeal

        theta = self.__ph.degToRad(theta)
        Co = self.__Clamb(nu,no)
        Ce = self.__Clamb(nu,ne)
        Po = self.__prop(nu, no, d)
        Pe = self.__prop(nu, ne, d)
        mPo = self.__mprop(nu, no, d)
        mPe = self.__mprop(nu, ne, d)
        M2invRet = 0.5*np.matrix([[Po*np.cos(theta),   Po*np.sin(theta),  -Po*np.sin(theta)/Co, Po*np.cos(theta)/Co],
                                  [-Pe*np.sin(theta),  Pe*np.cos(theta),  -Pe*np.cos(theta)/Ce, -Pe*np.sin(theta)/Ce],
                                  [mPo*np.cos(theta),  mPo*np.sin(theta), mPo*np.sin(theta)/Co, -mPo*np.cos(theta)/Co],
                                  [-mPe*np.sin(theta), mPe*np.cos(theta), mPe*np.cos(theta)/Ce, mPe*np.sin(theta)/Ce]])
        return M2invRet
    
    def __m(self, nu, theta, no1=None, ne1=None, no2=None, ne2=None, d=None):
        if no1 == None:
            no1 = self.__hwp.oN
        if ne1 == None:
            ne1 = self.__hwp.eN
        if no2 == None:
            no2 = self.__hwp.oN
        if ne2 == None:
            ne2 = self.__hwp.eN
        if d == None:
            d = self.__hwp.thickIdeal

        return self.__M1(nu,no1,ne1)*self.__M2inv(nu,theta,no2,ne2,d)

    def __M(self, nu, thetaArr=None, no=None, ne=None, d=None):
        if thetaArr == None:
            thetaArr = self.__hwp.interfaceAngles_3[1:]
        if no == None:
            no = self.__hwp.oN
        if ne == None:
            ne = self.__hwp.eN
        if d == None:
            d = self.__hwp.thickIdeal

        mRet = np.matrix([[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.]])
        for i in range(len(thetaArr)):
            mRet = mRet*self.__m(nu,thetaArr[i],no,ne,d)
        return mRet

    def __Atilde(self, nu, thetaArr=None, no=None, ne=None, d=None):
        if thetaArr == None:
            thetaArr = self.__hwp.interfaceAngles_3[1:]
        if no == None:
            no = self.__hwp.oN
        if ne == None:
            ne = self.__hwp.eN
        if d == None:
            d = self.__hwp.thickIdeal

        CC = self.__Clamb(nu, 1.)
        mm = self.__M(nu,thetaArr,no,ne,d)
        AtildeRet =  np.matrix([[mm.item(0,0)+CC*mm.item(0,3), mm.item(0,1)-CC*mm.item(0,2)],
                                [mm.item(1,0)+CC*mm.item(1,3), mm.item(1,1)-CC*mm.item(1,2)],
                                [mm.item(2,0)+CC*mm.item(2,3), mm.item(2,1)-CC*mm.item(2,2)],
                                [mm.item(3,0)+CC*mm.item(3,3), mm.item(3,1)-CC*mm.item(3,2)]])
        return AtildeRet

    def __Binv(self, nu, theta1, no=None, ne=None):
        if no == None:
            no = self.__hwp.oN
        if ne == None:
            ne = self.__hwp.eN

        CC = self.__Clamb(nu, 1.)
        theta1 = self.__ph.degToRad(theta1)
        BinvRet = 0.5*np.matrix([[np.cos(theta1),  np.sin(theta1), -np.sin(theta1)/CC, np.cos(theta1)/CC],
                                 [-np.sin(theta1), np.cos(theta1), -np.cos(theta1)/CC, -np.sin(theta1)/CC],
                                 [np.cos(theta1),  np.sin(theta1), np.sin(theta1)/CC,  -np.cos(theta1)/CC],
                                 [-np.sin(theta1), np.cos(theta1), np.cos(theta1)/CC,  np.sin(theta1)/CC]])
        return BinvRet

    def __a(self, nu, thetaArr=None, no=None, ne=None, d=None):
        if thetaArr == None:
            thetaArr = self.__hwp.interfaceAngles_3
        if no == None:
            no = self.__hwp.oN
        if ne == None:
            ne = self.__hwp.eN
        if d == None:
            d = self.__hwp.thickIdeal

        aRet = self.__Binv(nu, thetaArr[0], no, ne)*self.__Atilde(nu, thetaArr[1:], no, ne, d)
        return aRet

    def __r(self, nu, thetaArr=None, no=None, ne=None, d=None):
        if thetaArr == None:
            thetaArr = self.__hwp.interfaceAngles_3
        if no == None:
            no = self.__hwp.oN
        if ne == None:
            ne = self.__hwp.eN
        if d == None:
            d = self.__hwp.thickIdeal

        aa = self.__a(nu, thetaArr, no, ne, d)
        return np.matrix([[(aa.item(1,1)*aa.item(2,0)-aa.item(2,1)*aa.item(1,0))/(aa.item(1,1)*aa.item(0,0) - aa.item(0,1)*aa.item(1,0)), 
                           (aa.item(0,0)*aa.item(2,1)-aa.item(2,0)*aa.item(0,1))/(aa.item(1,1)*aa.item(0,0) - aa.item(0,1)*aa.item(1,0))],
                          [(aa.item(3,0)*aa.item(1,1)-aa.item(3,1)*aa.item(1,0))/(aa.item(1,1)*aa.item(0,0) - aa.item(0,1)*aa.item(1,0)), 
                           (aa.item(3,1)*aa.item(0,0)-aa.item(3,0)*aa.item(0,1))/(aa.item(1,1)*aa.item(0,0) - aa.item(0,1)*aa.item(1,0))]])

    def __t(self, nu, thetaArr=None, no=None, ne=None, d=None):
        if thetaArr == None:
            thetaArr = self.__hwp.interfaceAngles_3
        if no == None:
            no = self.__hwp.oN
        if ne == None:
            ne = self.__hwp.eN
        if d == None:
            d = self.__hwp.thickIdeal

        aa = self.__a(nu, thetaArr, no, ne, d)
        return np.matrix([[ aa.item(1,1)/(aa.item(1,1)*aa.item(0,0)-aa.item(0,1)*aa.item(1,0)), 
                            -aa.item(0,1)/(aa.item(1,1)*aa.item(0,0)-aa.item(0,1)*aa.item(1,0))],
                          [-aa.item(1,0)/(aa.item(1,1)*aa.item(0,0)-aa.item(0,1)*aa.item(1,0)), 
                            aa.item(0,0)/(aa.item(1,1)*aa.item(0,0)-aa.item(0,1)*aa.item(1,0))]])
    
    #***** Public functions *****
    #Transfer matrix for an ideal HWP
    def T(self, nu, thetaArr=None, no=None, ne=None, d=None):
        if thetaArr == None:
            thetaArr = self.__hwp.interfaceAngles_3
        if no == None:
            no = self.__hwp.oN
        if ne == None:
            ne = self.__hwp.eN
        if d == None:
            d = self.__hwp.thickIdeal

        tt = self.__t(nu, thetaArr, no=oN, ne=eN, d=singleThk)
        ts = np.conj(tt)
        return 0.5*np.matrix([[tt.item(0,0)*ts.item(0,0) + tt.item(1,0)*ts.item(1,0) + tt.item(0,1)*ts.item(0,1) + tt.item(1,1)*ts.item(1,1),
                               tt.item(0,0)*ts.item(0,0) + tt.item(1,0)*ts.item(1,0) - tt.item(0,1)*ts.item(0,1) - tt.item(1,1)*ts.item(1,1),
                               tt.item(0,0)*ts.item(0,1) + tt.item(1,0)*ts.item(1,1) + tt.item(0,1)*ts.item(0,0) + tt.item(1,1)*ts.item(1,0),
                               (1./1j)*(tt.item(0,0)*ts.item(0,1) + tt.item(1,0)*ts.item(1,1) - tt.item(0,1)*ts.item(0,0) - tt.item(1,1)*ts.item(1,0))],
                              
                              [tt.item(0,0)*ts.item(0,0) - tt.item(1,0)*ts.item(1,0) + tt.item(0,1)*ts.item(0,1) - tt.item(1,1)*ts.item(1,1),
                               tt.item(0,0)*ts.item(0,0) - tt.item(1,0)*ts.item(1,0) - tt.item(0,1)*ts.item(0,1) + tt.item(1,1)*ts.item(1,1),
                               tt.item(0,1)*ts.item(0,0) - tt.item(1,1)*ts.item(1,0) + tt.item(0,0)*ts.item(0,1) - tt.item(1,0)*ts.item(1,1),
                               (1./1j)*(tt.item(0,1)*ts.item(0,0) - tt.item(1,1)*ts.item(1,0) - tt.item(0,0)*ts.item(0,1) + tt.item(1,0)*ts.item(1,1))],
                              
                              [tt.item(0,0)*ts.item(1,0) + tt.item(1,0)*ts.item(0,0) + tt.item(0,1)*ts.item(1,1) + tt.item(1,1)*ts.item(0,1),
                               tt.item(0,0)*ts.item(1,0) + tt.item(1,0)*ts.item(0,0) - tt.item(0,1)*ts.item(1,1) - tt.item(1,1)*ts.item(0,1),
                               tt.item(0,0)*ts.item(1,1) + tt.item(1,0)*ts.item(0,1) + tt.item(0,1)*ts.item(1,0) + tt.item(1,1)*ts.item(0,0),
                               (1./1j)*(tt.item(0,0)*ts.item(1,1) + tt.item(1,0)*ts.item(0,1) - tt.item(0,1)*ts.item(1,0) - tt.item(1,1)*ts.item(0,0))],
                              
                              [1j*(tt.item(0,0)*ts.item(1,0) - tt.item(1,0)*ts.item(0,0) + tt.item(0,1)*ts.item(1,1) - tt.item(1,1)*ts.item(0,1)),
                               1j*(tt.item(0,0)*ts.item(1,0) - tt.item(1,0)*ts.item(0,0) - tt.item(0,1)*ts.item(1,1) + tt.item(1,1)*ts.item(0,1)),
                               1j*(tt.item(0,0)*ts.item(1,1) - tt.item(1,0)*ts.item(0,1) + tt.item(0,1)*ts.item(1,0) - tt.item(1,1)*ts.item(0,0)),
                               tt.item(0,0)*ts.item(1,1) - tt.item(1,0)*ts.item(0,1) - tt.item(0,1)*ts.item(1,0) + tt.item(1,1)*ts.item(0,0)]])
    
    #Output Stokes vector for an ideal HWP
    def Sout(self, alpha, nu, rho, P=1.):
        return self.T(nu,angles)*self.__ph.Stokes(P,alpha)
