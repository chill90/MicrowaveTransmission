#python Version 2.7.2
import numpy           as np
import scipy.constants as ct

#Class for calculating boundary conditions
class Hou:
    def __init__(self):
        self.__normInc = 0.
        return

    #Transmission through dielectric stack
    def calc(self, nArr, dArr, ltArr, freqArr, incAng=None):
        if incAng == None:
            incAng = self.__normInc
            
        #Number of interfaces to be analyzed
        if (len(nArr) == len(dArr)) and (len(nArr) == len(ltArr)):
            #Number of interfaces
            numInt = len(nArr) - 1
        else:
            raise Exception('Error in Hou.trans(): len(nArr), len(dArr), and len(ltArr) must be equal')

        #Calculate the angles of propogation in each layer
        theta = [incAng]
        for i in range(numInt):
            newAng = np.arcsin((nArr[i]/nArr[i+1])*np.sin(theta[i]))
            theta.append(newAng)
        theta = np.array(theta)
        #Calculate the reflection coefficients
        rs = np.array([(nArr[i]*np.cos(theta[i]) - nArr[i+1]*np.cos(theta[i+1]))/(nArr[i]*np.cos(theta[i]) + nArr[i+1]*np.cos(theta[i+1])) for i in range(numInt)])
        rp = np.array([(nArr[i]*(1./np.cos(theta[i])) - nArr[i+1]*(1./np.cos(theta[i+1])))/(nArr[i]*(1./np.cos(theta[i])) + nArr[i+1]*(1./np.cos(theta[i+1]))) for i in range(numInt)])
        #Calculate transmission coefficients
        ts = np.array([1.+rs[i] for i in range(0, numInt)])
        tp = np.array([1.+rp[i] for i in range(0, numInt)])
        
        #Calculate reflection and transmission for each frequency
        Ms = np.matrix([[1., rs[0]],[rs[0], 1.]])*(1./ts[0]) #Just the identity matrix for the s polarization
        Mp = np.matrix([[1., rp[0]],[rp[0], 1.]])*(1./tp[0]) #Just the identity matrix for the p polarization
        for i in range(1, numInt):
            #Calculate the exponential propogation factors
            pp = (2.*ct.pi*nArr[i]*dArr[i]*(1./np.cos(theta[i]))*np.array(freqArr)/ct.c)*( 0.5*ltArr[i] + 1.0j)
            pm = (2.*ct.pi*nArr[i]*dArr[i]*(1./np.cos(theta[i]))*np.array(freqArr)/ct.c)*(-0.5*ltArr[i] - 1.0j)
            #Calculate new propogation matrices
            Ms = Ms*np.matrix([[np.exp(pp), 0],[0, np.exp(pm)]])*np.matrix([[1., rs[i]],[rs[i], 1.]])*(1./ts[i])
            Mp = Mp*np.matrix([[np.exp(pp), 0],[0, np.exp(pm)]])*np.matrix([[1., rp[i]],[rp[i], 1.]])*(1./tp[i])
            
        #Frequencies
        Freq = freqArr
        
        #Calculate transmitted power
        Tran_s = np.array(list(map(lambda x: abs(x)**2, 1./(Ms.item(0,0)))))
        Tran_p = np.array(list(map(lambda x: abs(x)**2, 1./(Mp.item(0,0)))))

        #Calculate reflected power
        Refl_s = np.array(list(map(lambda x: abs(x)**2, Ms.item(1,0)/Ms.item(0,0))))
        Refl_p = np.array(list(map(lambda x: abs(x)**2, Mp.item(1,0)/Mp.item(0,0))))

        #Calculate absorbed power
        Abso_s = 1. - Tran_s - Refl_s
        Abso_p = 1. - Tran_p - Refl_p
        
        return Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s
