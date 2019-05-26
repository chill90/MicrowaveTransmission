import numpy as np
import          os

#Class for importing measured data from various measurement setups
class Measurement:
    def __init__(self, dataFile=None):
        #Store passed parameters
        self.dataFile = dataFile
        self.fhandle  = self.dataFile.split('.')[0].split('/')[-1]
        #Allowed measurement apparatuses
        self.app = {'michRefl': 'UMICH_REFLECTOMETER',
                    'dick': 'DICK_COHERENTSOURCE'}
        
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
            self.freq = self.freq
            self.ptrans  = None;    self.ptranserr = None
            self.strans  = None;    self.stranserr = None  
            self.prefl   = refl/2.; self.preflerr = None
            self.srefl   = refl/2.; self.sreflerr = None
            self.pabso   = None;    self.pabsoerr = None
            self.sabso   = None;    self.sabsoerr = None
        if self.app['dick'] in dataFile.upper():
            freq, trans, transerr = np.loadtxt(dataFile, unpack=True, comments='#', usecols=[0,4,5])
            self.freq = freq;
            self.ptrans  = trans; self.ptranserr = transerr/np.sqrt(2.)
            self.strans  = trans; self.stranserr = transerr/np.sqrt(2.)
            self.prefl   = None;  self.preflerr = None
            self.srefl   = None;  self.sreflerr = None
            self.pabso   = None;  self.pabsoerr = None
            self.sabso   = None;  self.sabsoerr = None
        else:
            raise Exception("MICROWAVE TRANSMISSION ERROR: Unable to find valid instrument descriptor in data file name '%s'" % (dataFile))
        self.outputs = (self.freq, 
                        self.ptrans, self.ptranserr,
                        self.strans, self.stranserr,
                        self.prefl,  self.preflerr, 
                        self.srefl,  self.sreflerr,
                        self.pabso,  self.pabsoerr,
                        self.sabso,  self.sabsoerr)
