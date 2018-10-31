#Using python 2.7.2
import numpy             as np
import matplotlib.pyplot as plt
import                      os

#Custom classes
import src.measurement   as ms

class Plot:
    def __init__(self, saveLoc=None, fhandle=None, flo=None, fhi=None):
        if fhandle is None:
            self.fname = ''
        else:
            self.fname = '_%s' % (fhandle)
        if saveLoc is None:
            self.saveLoc = '..'+os.sep+'Plots'+os.sep
        else:
            self.saveLoc = saveLoc+os.sep
        #Measurement handles
        self.meas        = ms.Measurement()
        self.measHandles = self.meas.app.values()
            
        #Set initial plotting parameters
        plt.rc('font', family='serif')
        plt.rc('font', size=28)
        self.lw = 4
        self.fignum = 0
        self.lfsz = 18
        self.ms   = 6
        self.figsize = (16,12)
        
    #Method for plotting output data from transmission calculation
    def plotTrans(self, outputs, labels):
        if not len(outputs) == len(labels):
            raise Exception("MICROWAVE TRANSMISSION EXCEPTION: In plot.plotRefl(), number of outputs and number of labels have different lengths: %d, %d" % (len(outputs), len(labels)))
        plt.figure(self.fignum, figsize=self.figsize)
        for i in range(len(outputs)):
            output = outputs[i]
            label  = labels[i]
            #Plot central values
            if label.upper() in self.measHandles:
                self.__datPlot(output[0], np.mean([output[1], output[3]], axis=0), np.sqrt(np.power(output[2],2)+np.power(output[4],2)), label=label)
            else:
                self.__simPlot(output[0], np.mean([output[1], output[3]], axis=0), np.sqrt(np.power(output[2],2)+np.power(output[4],2)), label=label)
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Transmission')
        #plt.ylim([0,1])
        self.__legend()
        plt.savefig('%sTransmission%s.jpg' % (self.saveLoc, self.fname))
        self.fignum += 1

#   def plotTrans(self, output):
#       plt.figure(self.fignum, figsize=self.figsize)
#       #Plot central values
#       self.__plot(output[0], np.mean([output[1], output[3]], axis=0), np.sqrt(np.power(output[2],2)+np.power(output[4],2)), color='r', label='Simulation')
#       plt.xlabel('Frequency [GHz]')
#       plt.ylabel('Transmission')
#       plt.ylim([0,1])
#       self.__legend()
#       plt.savefig('%sTransmission%s.jpg' % (self.saveLoc, self.fname))
#       self.fignum += 1
#
    def plotRefl(self, outputs, labels):
        if not len(outputs) == len(labels):
            raise Exception("MICROWAVE TRANSMISSION EXCEPTION: In plot.plotRefl(), number of outputs and number of labels have different lengths: %d, %d" % (len(outputs), len(labels)))
        plt.figure(self.fignum, figsize=self.figsize)
        for i in range(len(outputs)):
            output = outputs[i]
            label  = labels[i]
            #Plot central values
            if label.upper() in self.measHandles:
                self.__datPlot(output[0], np.mean([output[5], output[7]], axis=0), np.sqrt(np.power(output[6],2)+np.power(output[8],2)), label=label)
            else:
                self.__simPlot(output[0], np.mean([output[5], output[7]], axis=0), np.sqrt(np.power(output[6],2)+np.power(output[8],2)), label=label)
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Reflection')
        plt.ylim([0,1])
        self.__legend()
        plt.savefig('%sReflection%s.jpg' % (self.saveLoc, self.fname))
        self.fignum += 1

    def plotAbso(self, output):
        plt.figure(self.fignum, figsize=self.figsize)
        #Plot central values
        self.__plot(output[0], np.mean([output[9], output[11]], axis=0), np.sqrt(np.power(output[10],2)+np.power(output[12],2)), color='k', label='Simulation')
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Absorption')
        plt.ylim([0,1])
        self.__legend()
        plt.savefig('%sAbsorption%s.jpg' % (self.saveLoc, self.fname))
        self.fignum += 1
        
    def plotAll(self, output):
        self.plotTrans(output)
        self.plotRefl(output)
        self.plotAbso(output)
        return

    #Generic funciton for plotting
    def __simPlot(self, freq, cent, err, label):
        upper = np.array(cent) + np.array(err)
        lower = np.array(cent) - np.array(err)
        p = plt.plot(freq, cent, linewidth=self.lw, label=label)
        plt.fill_between(freq, lower, upper, color=p[0].get_color(), linewidth=0, alpha=0.25)
        return 1

    #Function for overplotting data from various apparatuses
    def __datPlot(self, freq, cent, err, label):
        plt.errorbar(freq, cent, yerr=err, marker='o', markersize=self.ms, linewidth=self.lw, label=label)
        return 1

    #Generic function for generting legend
    def __legend(self):
        plt.legend(loc='best', fontsize=self.lfsz)
        return 1
    
