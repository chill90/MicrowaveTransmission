#Using python 2.7.2
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
                          "Sim Method":                              str(values[params.index("Sim Method")]                                                )}   
        
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
        return (self.freqs*un.Hz_to_GHz,
                np.mean(tran_p, axis=0), np.std(tran_p, axis=0),
                np.mean(tran_s, axis=0), np.std(tran_s, axis=0),
                np.mean(refl_p, axis=0), np.std(refl_p, axis=0), 
                np.mean(refl_s, axis=0), np.std(refl_s, axis=0), 
                np.mean(abso_p, axis=0), np.std(abso_p, axis=0),
                np.mean(abso_s, axis=0), np.std(abso_s, axis=0))
        
