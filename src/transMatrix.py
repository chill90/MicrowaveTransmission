#!/usr/local/bin/python

#python Version 2.7.2
import numpy as np
import glob
import pickle as pkl
import scipy.integrate as itg
import copy

#Class to calculate transmission through arbitrary stack
#Method taken from Tom Essinger-Hileman's paper from 2014
#"Transfer matrix for treating stratified media including birefringent crystals"
class Tom:
    def __init__(self):
        #Private objects
        self.__ph = Physics()
        self.__hwp = HWP()

        #Private variables
        self.__defHWPAng = 0.
        self.__defIncAng = 0.
        self.__defPolAng = 0.

        return 

    #***** Private functions *****
    #Invert matrix
    def __inv(self, mat):
        return np.linalg.inv(mat)
    #Degree to radian
    def __dr(self, theta):
        return self.__ph.degToRad(theta)
    #Radian to degree
    def __rd(self, theta):
        return self.__ph.radToDeg(theta)
    #Cosine
    def __cos(self, theta):
        return np.cos(self.__dr(theta))
    #Sine
    def __sin(self, theta):
        return np.sin(self.__dr(theta))
    #Exponent
    def __pow(self, arg, index):
        return np.power(arg, index)
    #Square
    def __sq(self, arg):
        return np.power(arg, 2.)

    #Transfer function for stratified medium
    def TransMat(self, nu, noArr, neArr, oLTArr, eLTArr, tArr, chiArr, rho, incAngle, incN):
        #First medium is air
        theta1 = incAngle
        n1 = incN
    
        #Create transfer matrix
        Tr = np.identity(4)
        #Return the identity matrix if no layers are to be calculated
        if len(noArr) == 0:
            return Tr
        #Loop over the HWP layers
        for i in range(len(chiArr)):
            #Physical values for this layer
            no = noArr[i]
            ne = neArr[i]
            oLT = oLTArr[i]
            eLT = eLTArr[i]
            t = tArr[i]
            chi = chiArr[i]+rho #Plate orientation w.r.t x-axis
            
            #Vacuum wavevector
            k0 = (2*np.pi)/(self.__ph.c/nu)
            
            nP = no
            nPP = ne*np.sqrt(1+(self.__pow(ne,-2)-self.__pow(no,-2))*self.__sq(n1)*self.__sq(self.__sin(theta1))*self.__sq(self.__cos(chi)))
        
            R = lambda chi: np.matrix([[self.__cos(chi),-self.__sin(chi),0],[self.__sin(chi),self.__cos(chi),0],[0,0,1]])
            epsilonP = R(chi)*np.matrix([[self.__sq(ne),0,0],[0,self.__sq(no),0],[0,0,self.__sq(no)]])*R(-chi)
        
            LTp = oLT
            nPT = nP*np.sqrt(1-1j*LTp)
            
            if eLT == oLT:
                LTpp = oLT
            else:
                LTpp = oLT + ((eLT - oLT)/(ne - no))*(nPP - no)
            nPPT = nPP*np.sqrt(1-1j*LTpp)
            
            thetaP = self.__rd(np.arcsin(n1*self.__sin(theta1)/nP))
            thetaPP = self.__rd(np.arcsin(n1*self.__sin(theta1)/nPP))
            
            deltaP = nPT*t*self.__cos(thetaP)
            deltaPP = nPPT*t*self.__cos(thetaPP)
            
            DPt1 = self.__pow(self.__sq(self.__cos(thetaP))+self.__sq(self.__sin(thetaP))*self.__sq(self.__sin(chi)),-0.5)*np.array([-self.__sin(chi)*self.__cos(thetaP),self.__cos(chi)*self.__cos(thetaP),self.__sin(chi)*self.__sin(thetaP)])
            DPPt1 = self.__pow(self.__sq(self.__cos(chi))*self.__sq(self.__cos(thetaP))+self.__sq(self.__sin(chi))*self.__sq(self.__cos(thetaP-thetaPP)),-0.5)*np.array([self.__cos(chi)*self.__cos(thetaP)*self.__cos(thetaPP),self.__sin(chi)*(self.__sin(thetaP)*self.__sin(thetaPP)+self.__cos(thetaP)*self.__cos(thetaPP)),-self.__cos(chi)*self.__cos(thetaP)*self.__sin(thetaPP)])
            HPt1 = self.__pow(self.__sq(self.__cos(thetaP))*self.__sq(self.__cos(chi))+self.__sq(self.__sin(chi)),-0.5)*np.array([-self.__sq(self.__cos(thetaP))*self.__cos(chi),-self.__sin(chi),self.__cos(thetaP)*self.__sin(thetaP)*self.__cos(chi)])
            HPPt1 = self.__pow(self.__sq(self.__cos(thetaP-thetaPP))*self.__sq(self.__sin(chi))+self.__sq(self.__cos(thetaP))*self.__sq(self.__cos(chi)),-0.5)*np.array([-self.__cos(thetaP-thetaPP)*self.__cos(thetaPP)*self.__sin(chi),self.__cos(thetaP)*self.__cos(chi),self.__cos(thetaP-thetaPP)*self.__sin(thetaPP)*self.__sin(chi)])
        
            Phi1 = np.matrix([[DPt1.item(0),DPPt1.item(0),DPt1.item(0),DPPt1.item(0)],[(1./nP)*HPt1.item(1),(1./nPP)*HPPt1.item(1),(-1./nP)*HPt1.item(1),(-1./nPP)*HPPt1.item(1)],[DPt1.item(1),DPPt1.item(1),DPt1.item(1),DPPt1.item(1)],[(-1./nP)*HPt1.item(0),(-1./nPP)*HPPt1.item(0),(1./nP)*HPt1.item(0),(1./nPP)*HPPt1.item(0)]])
            Psi1 = np.matrix([[self.__inv(epsilonP).item(0,0),0,self.__inv(epsilonP).item(0,1),0],[0,1.,0,0],[self.__inv(epsilonP).item(1,0),0,self.__inv(epsilonP).item(1,1),0],[0,0,0,1.]])
            DeltaP = 1j*k0*deltaP
            DeltaPP = 1j*k0*deltaPP
            P = np.matrix([[np.exp(-DeltaP),0,0,0],[0,np.exp(-DeltaPP),0,0],[0,0,np.exp(DeltaP),0],[0,0,0,np.exp(DeltaPP)]])
        
            #Store the transfer matrix fragment into an array
            Tr = Psi1*Phi1*self.__inv(P)*self.__inv(Phi1)*self.__inv(Psi1)*Tr
        
        #Return the transfer matrix
        return Tr

    #Transmission through dielectric layers
    def trans(self, nuArr, noArr=None, neArr=None, oLTArr=None, eLTArr=None, tArr=None, chiArr=None, rho=None, incAng=None, polAng=None):
        if noArr == None:
            noArr = self.__hwp.noArr
        if neArr == None:
            neArr = self.__hwp.neArr
        if oLTArr == None:
            oLTArr = self.__hwp.oltArr
        if eLTArr == None:
            eLTArr = self.__hwp.eltArr
        if tArr == None:
            tArr = self.__hwp.dArr
        if chiArr == None:
            chiArr = self.__hwp.plateAngles
        if rho == None:
            rho = self.__defHWPAng
        if incAng == None:
            incAng = self.__defIncAng
        if polAng == None:
            polAng = self.__defPolAng
        
        #Check array integrity
        if not (len(chiArr) == len(noArr) == len(neArr) == len(tArr) == len(oLTArr) == len(eLTArr)):
            raise NameError("Length of 'anisotropicTransmission' arrays need to be the same")
        if not len(nuArr):
            raise NameError("Length 'nuArr' needs to be non-zero")

        #Cannot calculate birefringence in first layer, so take average index
        n1 = (noArr[0] + neArr[0])/2.
        theta1 = incAng
        
        #Create transfer matrix for this stack
        T = lambda nu: self.TransMat(nu, noArr[1:-1], neArr[1:-1], oLTArr[1:-1], eLTArr[1:-1], tArr[1:-1], chiArr[1:-1], rho, incAng, n1)
        
        #Cannot calculate birefringence in final layer, so take average index
        n3 = (noArr[-1] + neArr[-1])/2.
        theta3 = self.__rd(np.arcsin((n1/n3)*self.__sin(theta1)))
        
        #Transmission coefficients
        alphaT = lambda nu: (T(nu).item(0,0)*self.__cos(theta3)+T(nu).item(0,1)*n3)/self.__cos(theta1)
        betaT = lambda nu: (T(nu).item(0,2)+T(nu).item(0,3)*n3*self.__cos(theta3))/self.__cos(theta1)
        gammaT = lambda nu: (T(nu).item(1,0)*self.__cos(theta3)+T(nu).item(1,1)*n3)/n1
        deltaT = lambda nu: (T(nu).item(1,2)+T(nu).item(1,3)*n3*self.__cos(theta3))/n1
        etaT = lambda nu: (T(nu).item(2,0)*self.__cos(theta3)+T(nu).item(2,1)*n3)
        kappaT = lambda nu: (T(nu).item(2,2)+T(nu).item(2,3)*n3*self.__cos(theta3))
        rhoT = lambda nu: (T(nu).item(3,0)*self.__cos(theta3)+T(nu).item(3,1)*n3)/(n1*self.__cos(theta1))
        sigmaT = lambda nu: (T(nu).item(3,2)+T(nu).item(3,3)*n3*self.__cos(theta3))/(n1*self.__cos(theta1))
        GammaT = lambda nu: self.__pow((alphaT(nu)+gammaT(nu))*(kappaT(nu)+sigmaT(nu))-(betaT(nu)+deltaT(nu))*(etaT(nu)+rhoT(nu)),-1.)
        
        #Jones Transmission and Reflection Matrices
        JT = lambda nu: 2.*GammaT(nu)*np.matrix(
            [[(kappaT(nu)+sigmaT(nu)),
              (-betaT(nu)-deltaT(nu))],
             [(-etaT(nu)-rhoT(nu)),
              (alphaT(nu)+gammaT(nu))]]
            )
        JR = lambda nu: GammaT(nu)*np.matrix(
            [[(gammaT(nu)-alphaT(nu))*(kappaT(nu)+sigmaT(nu))-(deltaT(nu)-betaT(nu))*(etaT(nu)+rhoT(nu)),
              2.*(alphaT(nu)*deltaT(nu)-gammaT(nu)*betaT(nu))],
             [2.*(etaT(nu)*sigmaT(nu)-rhoT(nu)*kappaT(nu)),
              (alphaT(nu)+gammaT(nu))*(kappaT(nu)-sigmaT(nu))-(betaT(nu)+deltaT(nu))*(etaT(nu)-rhoT(nu))]]
            )
        
        #Incident Electric Field Amplitude
        Ei = np.array([abs(self.__cos(polAng)), abs(self.__sin(polAng))])
        #Incident Power
        Pi = np.array([self.__sq(Ei[0]), self.__sq(Ei[1])])
        
        #Transmitted Electric Field Amplitude
        EPT = lambda nu: JT(nu).item(0,0)*Ei[0] + JT(nu).item(0,1)*Ei[1]
        EST = lambda nu: JT(nu).item(1,0)*Ei[0] + JT(nu).item(1,1)*Ei[1]
        
        #Transmitted Power
        PPT = lambda nu: EPT(nu)*np.conj(EPT(nu))
        PST = lambda nu: EST(nu)*np.conj(EST(nu))

        #Reflected Electric Field Amplitude
        EPR = lambda nu: JR(nu).item(0,0)*Ei[0] + JR(nu).item(0,1)*Ei[1]
        ESR = lambda nu: JR(nu).item(1,0)*Ei[0] + JR(nu).item(1,1)*Ei[1]
        
        #Reflected Power
        PPR = lambda nu: EPR(nu)*np.conj(EPR(nu))
        PSR = lambda nu: ESR(nu)*np.conj(ESR(nu))
        
        #Loop over the frequencies and calculate transmission for each
        pT = []
        sT = []
        pR = []
        sR = []
        for nu in nuArr:
            #print "Analyzing %.1f GHz" % (nu*1.e-9)
            pT.append(PPT(nu))
            sT.append(PST(nu))
            pR.append(PPR(nu))
            sR.append(PSR(nu))

        #Absorbed power
        pA = 1. - np.array(pT) - np.array(pR)
        sA = 1. - np.array(sT) - np.array(sR)
            
        #Return the power transmission, reflection, and absorption arrays
        return np.array(pT), np.array(sT), np.array(pR), np.array(sR), np.array(pA), np.array(sA)
