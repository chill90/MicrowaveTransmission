import numpy as np
import          os

#Class for importing measured data from various measurement setups
class Measurement:
    def __init__(self, dataFile=None):
        #Store passed parameters
        self.dataFile = dataFile
        #Allowed measurement apparatuses
        self.app = {'michRefl': 'UMICH_REFLECTOMETER'}
        
    # ***** Public Methods *****
    def loadData(self, dataFile=None):
        if dataFile is None:
            dataFile = self.dataFile
        if dataFile is None:
            raise Exception("MICROWAVE TRANSMISSION ERROR: No data file to be processes for Measurement class")
        if not os.path.isfile(dataFile):
            raise Exception("MICROWAVE TRANSMISSION ERROR: Unable to locate data file '%s'" % (dataFile))
        if self.app['michRefl'] in dataFile.upper():
            freq, refl = np.loadtxt(dataFile, unpack=True, comments='#')
            freqerr = None
            ptrans  = None; ptranserr = None
            strans  = None; stranserr = None  
            prefl   = refl/2.; preflerr = None
            srefl   = refl/2.; sreflerr = None
            pabso   = None;    pabsoerr = None
            sabso   = None;    sabsoerr = None
        else:
            raise Exception("MICROWAVE TRANSMISSION ERROR: Unable to find valid instrument descriptor in data file name '%s'" % (dataFile))
        return (freq, 
                ptrans, ptranserr,
                strans, stranserr,
                prefl,  preflerr, 
                srefl,  sreflerr,
                pabso,  pabsoerr,
                sabso,  sabsoerr)
