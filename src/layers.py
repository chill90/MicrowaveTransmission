#Using python 2.7.2
import numpy       as np
import collections as cl

#Custom classes
import src.parameter  as pm
import src.unit       as un

class Layers:
    def __init__(self, file):
        #Store input variables
        self.file = file
        
        #Parse file for layer parameters
        layers, thicks, indexes, lossTans = np.loadtxt(self.file, comments='+-', skiprows=2, delimiter='|', dtype=np.str, unpack=True)[1:-1]
        #Build dictionaries of dielectric layers
        self.layers = cl.OrderedDict({})
        for i in range(len(layers)):
            self.layers[layers[i]] = {"Thickness": pm.Parameter("Thickness", thicks[i],   unit=un.mm_to_m, min=0.0, max=np.inf),
                                      "Index":     pm.Parameter("Index",     indexes[i],                   min=0.0, max=np.inf),
                                      "LossTan":   pm.Parameter("LossTan",   lossTans[i], unit=1.e-04,     min=0.0, max=np.inf)}
        
    #Method to sample layers
    def sample(self):
        thicks = []; indexes = []; lossTans = []
        for k in self.layers:
            thicks.append(  self.layers[k]["Thickness"].sample())
            indexes.append( self.layers[k]["Index"    ].sample())
            lossTans.append(self.layers[k]["LossTan"  ].sample())
        return thicks, indexes, lossTans

    #Method to get average values for layers
    def getAvg(self):
        thicks = []; indexes = []; lossTans = []
        for k in self.layers:
            thicks.append(  self.layers[k]["Thickness"].getAvg())
            indexes.append( self.layers[k]["Index"    ].getAvg())
            lossTans.append(self.layers[k]["LossTan"  ].getAvg())
        return thicks, indexes, lossTans
    
