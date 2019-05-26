#Using> python 2.7.2
import numpy       as np
import sys         as sy
import                os

#Custom classes
import src.layers     as ly
import src.parameter  as pm
import src.unit       as un
import src.hou        as ho

class Simulate:
    def __init__(self, layerFile=None, simFile=None):
        #Generate layers instance
        if layerFile is None:
            self.layerFile = os.path.abspath(__file__)+'..'+os.sep+'config'+os.sep+'layers'+os.sep+'layers.txt'
        else:
            self.layerFile = layerFile
        self.fhandle = self.layerFile.split('.')[0].split('/')[-1]
        self.layers = ly.Layers(self.layerFile)

        #Save input parameters
        if simFile is None:
            self.simFile = os.path.abspath(__file__)+'..'+os.sep+'config'+os.sep+'simulation'+os.sep+'simInputs.txt'
        else:
            self.simFile = simFile
        #Load simulation input parameters
        params, values = np.loadtxt(self.simFile, unpack=True, usecols=[1,3], comments='+-', delimiter='|', dtype=np.str)
        params = [param.strip() for param in params]; values = [value.strip() for value in values]
        self.simInputs = {"Inc Angle":     pm.Parameter("Inc Angle",     values[params.index("Inc Angle")],     min=-90.,    max=90.,    unit=un.deg_to_rad),
                          "Inc Pol Angle": pm.Parameter("Inc Pol Angle", values[params.index("Inc Pol Angle")], min=-np.inf, max=np.inf, unit=un.deg_to_rad),
                          "Low Freq":      pm.Parameter("Low Freq",      values[params.index("Low Freq")],      min=0.0,     max=np.inf, unit=un.GHz_to_Hz ),
                          "High Freq":     pm.Parameter("High Freq",     values[params.index("High Freq")],     min=0.0,     max=np.inf, unit=un.GHz_to_Hz ),
                          "Freq Step":     pm.Parameter("Freq Step",     values[params.index("Freq Step")],     min=1.e-4,   max=1.e12,  unit=un.GHz_to_Hz ),
                          "Num Trials":                              int(values[params.index("Num Trials")]                                                ),
                          "Sim Method":                              str(values[params.index("Sim Method")]                                                ),
                          "Band Centers":  pm.Parameter("Band Centers", values[params.index("Band Centers")],  min=0.0,     max=np.inf, unit=un.GHz_to_Hz ),   
                          "Bandwidths":    pm.Parameter("Bandwidths",   values[params.index("Bandwidths")],    min=0.0,     max=2.0                       )}
        
        #Set frequency array
        self.freqs = np.arange(self.simInputs["Low Freq"].getAvg(), self.simInputs["High Freq"].getAvg()+self.simInputs["Freq Step"].getAvg(), self.simInputs["Freq Step"].getAvg())
    
    #Run simulation
    def calc(self):
        if self.simInputs["Sim Method"].upper() == 'HOU':
            return self.houCalc()
        else:
            print ("Error processing calculation in simulate.py")
            return 0
        
    #Simulate using Hou code
    def houCalc(self):
        #Save transmission, reflection, and absorption
        tran_p = []; tran_s = []
        refl_p = []; refl_s = []
        abso_p = []; abso_s = []
        #Instantiate Hou object
        hou = ho.Hou()
        #Simulate over the specified number of iterations
        for i in range(self.simInputs["Num Trials"]):
            if not i:
                #Sample the layers and incident angle
                thicks, indexes, lossTans = self.layers.getAvg()
                incAngle = self.simInputs["Inc Angle"].getAvg()
            else:
                #Sample the layers and incident angle
                thicks, indexes, lossTans = self.layers.sample()
                incAngle = self.simInputs["Inc Angle"].sample()                
            #Calculate the transmission
            freq, t_p, t_s, r_p, r_s, a_p, a_s = hou.calc(indexes, thicks, lossTans, self.freqs, incAngle)
            tran_p.append(t_p); tran_s.append(t_s)
            refl_p.append(r_p); refl_s.append(r_s)
            abso_p.append(a_p); abso_s.append(a_s)
        #Return all means and standard deviations
        self.freqs = self.freqs
        self.tran_p = self.__median(tran_p)
        self.tran_p_5 = self.__perc(tran_p, 5.)
        self.tran_p_95 = self.__perc(tran_p, 95.)
        self.tran_s = self.__median(tran_s)
        self.tran_s_5 = self.__perc(tran_s, 5.)
        self.tran_s_95 = self.__perc(tran_s, 95.)
        self.refl_p = self.__median(refl_p)
        self.refl_p_5 = self.__perc(refl_p, 5.)
        self.refl_p_95 = self.__perc(refl_p, 95.)
        self.refl_s = self.__median(refl_s)
        self.refl_s_5 = self.__perc(refl_s, 5.)
        self.refl_s_95 = self.__perc(refl_s, 95.)
        self.abso_p = self.__median(abso_p)
        self.abso_p_5 = self.__perc(abso_p, 5.)
        self.abso_p_95 = self.__perc(abso_p, 95.)
        self.abso_s = self.__median(abso_s)
        self.abso_s_5 = self.__perc(abso_s, 5.)
        self.abso_s_95 = self.__perc(abso_s, 95.)
        
        self.outputs = (self.freqs*un.Hz_to_GHz,
                        self.tran_p, self.tran_p_5, self.tran_p_95,
                        self.tran_s, self.tran_s_5, self.tran_s_95,
                        self.refl_p, self.refl_p_5, self.refl_p_95,
                        self.refl_s, self.refl_s_5, self.refl_s_95,
                        self.abso_p, self.abso_p_5, self.abso_p_95,
                        self.abso_s, self.abso_s_5, self.abso_s_95)
        return True
    
    #Function to find the percentile
    def __median(self, arr):
        if self.simInputs["Num Trials"] == 1:
            return arr[0]
        else:
            bas = self.__ba(arr)
            ba_med = np.median(bas)
            ind = np.argmin(abs(bas - ba_med))
            #ind = bas.index(ba_med)
            return arr[ind]

    #Function to find the percentile
    def __perc(self, arr, perc):
        if self.simInputs["Num Trials"] == 1:
            return np.array(len(arr[0])*[0.])
        else:
            bas = self.__ba(arr)
            ba_per = np.percentile(bas, perc)
            ind = np.argmin(abs(bas - ba_per))
            return arr[ind]

    #Function to find the band-averaged reflection
    def __ba(self, arr):
        self.bandCenters = np.array([89.5, 147.8])*1.e9
        self.bandwidths  = np.array([0.324,0.260])
        band_lo = self.bandCenters*(1. - 0.5*self.bandwidths)
        band_hi = self.bandCenters*(1. + 0.5*self.bandwidths)
        #Find average of each value within each band
        bas = []
        for i in range(len(arr)):
            mask = np.array([False]*len(self.freqs))
            for j in range(len(band_lo)):
                mask = np.logical_or((self.freqs < band_hi[j])*(self.freqs > band_lo[j]), mask)
            ba = np.mean(arr[i][mask])
            bas.append(ba)
        return bas
