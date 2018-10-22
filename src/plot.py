#Using python 2.7.2
import numpy             as np
import matplotlib.pyplot as plt
import                      os

class Plot:
    def __init__(self, saveLoc=None, fhandle=None):
        if fhandle is None:
            self.fname = ''
        else:
            self.fname = '_%s' % (fhandle)
        if saveLoc is None:
            self.saveLoc = '..'+os.sep+'Plots'+os.sep
        else:
            self.saveLoc = saveLoc+os.sep
        #Set initial plotting parameters
        plt.rc('font', family='serif')
        plt.rc('font', size=28)
        self.lw = 4
        self.fignum = 0
        self.lfsz = 18
        self.figsize = (16,12)
        
    #Method for plotting output data from transmission calculation
    def plotTrans(self, output):
        plt.figure(self.fignum, figsize=self.figsize)
        #Plot central values
        self.__plot(output[0], np.mean([output[1], output[3]], axis=0), np.sqrt(np.power(output[2],2)+np.power(output[4],2)), color='r', label='Transmission')
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Transmission')
        plt.savefig('%sTransmission%s.jpg' % (self.saveLoc, self.fname))
        self.fignum += 1

    def plotRefl(self, output):
        plt.figure(self.fignum, figsize=self.figsize)
        #Plot central values
        self.__plot(output[0], np.mean([output[5], output[7]], axis=0), np.sqrt(np.power(output[6],2)+np.power(output[8],2)), color='b', label='Reflection')
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Reflection')
        plt.savefig('%sReflection%s.jpg' % (self.saveLoc, self.fname))
        self.fignum += 1

    def plotAbso(self, output):
        plt.figure(self.fignum, figsize=self.figsize)
        #Plot central values
        self.__plot(output[0], np.mean([output[9], output[11]], axis=0), np.sqrt(np.power(output[10],2)+np.power(output[12],2)), color='k', label='Absorption')
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Absorption')
        plt.savefig('%sAbsorption%s.jpg' % (self.saveLoc, self.fname))
        self.fignum += 1
        
    def plotAll(self, output):
        self.plotTrans(output)
        self.plotRefl(output)
        self.plotAbso(output)
        return

    #Generic funciton for plotting
    def __plot(self, freq, cent, err, label, color):
        upper = np.array(cent) + np.array(err)
        lower = np.array(cent) - np.array(err)
        plt.plot(freq, cent, color=color, linewidth=self.lw, label=label)
        plt.fill_between(freq, lower, upper, color=color, linewidth=0, alpha=0.25)
        return 1
    
    #Generic function for generting legend
    def __legend(self):
        plt.legend(loc=best, fontsize=self.lfsz)
        return 1
    
