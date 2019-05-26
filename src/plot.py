#Using python 2.7.2
import numpy             as np
import matplotlib.pyplot as plt
import                      os

#Custom classes
import src.measurement   as ms
import src.unit          as un

class Plot:
    def __init__(self, sims, dats, saveLoc=None, opBands=True):
        self.sims = sims
        self.dats = dats
        self.fhandle = '_'.join([sim.fhandle for sim in self.sims] + [dat.fhandle for dat in self.dats])
        if saveLoc is None:
            self.saveLoc = '..'+os.sep+'Plots'+os.sep
        else:
            self.saveLoc = saveLoc+os.sep
        self.opBands = opBands
        #Measurement handles
        #self.meas        = ms.Measurement()
        #self.measHandles = self.meas.app.values()
            
        #Set initial plotting parameters
        plt.rc('font', family='serif')
        plt.rc('font', size=28)
        self.lw = 4
        self.fignum = 0
        self.lfsz = 18
        self.ms   = 6
        self.figsize = (16,12)

        #Bands to overplot
        if opBands:
            self.bandCenters = np.array(self.sims[0].simInputs["Band Centers"].avg)*un.Hz_to_GHz
            self.bandwidths = np.array(self.sims[0].simInputs["Bandwidths"].avg)
            self.bandLo = self.bandCenters*(1. - 0.5*self.bandwidths)
            self.bandHi = self.bandCenters*(1. + 0.5*self.bandwidths)
        
    #Method for plotting output data from transmission calculation
    def plotTrans(self):
        plt.figure(self.fignum, figsize=self.figsize)
        for sim in self.sims:
            self.__simPlot(sim.outputs[0], 
                           np.mean([sim.outputs[1], sim.outputs[4]], axis=0), 
                           np.mean([sim.outputs[2], sim.outputs[5]], axis=0), 
                           np.mean([sim.outputs[3], sim.outputs[6]], axis=0), 
                           label=sim.fhandle)
            #Calculate the band-averaged transmission
            for i in range(len(self.bandLo)):
                mask = (sim.outputs[0] > self.bandLo[i])*(sim.outputs[0] < self.bandHi[i])
                mean_ba = np.trapz(np.mean([sim.outputs[1][mask], sim.outputs[4][mask]], axis=0), sim.outputs[0][mask])/(sim.outputs[0][mask][-1] - sim.outputs[0][mask][0])
                lo_ba  = np.trapz(np.mean([sim.outputs[2][mask], sim.outputs[5][mask]], axis=0), sim.outputs[0][mask])/(sim.outputs[0][mask][-1] - sim.outputs[0][mask][0])
                hi_ba  = np.trapz(np.mean([sim.outputs[3][mask], sim.outputs[6][mask]], axis=0), sim.outputs[0][mask])/(sim.outputs[0][mask][-1] - sim.outputs[0][mask][0])
                print ('%.1f GHz Band Transmission for %s = %.3f + %.3f / - %.3f' % (self.bandCenters[i], sim.fhandle, mean_ba, hi_ba-mean_ba, mean_ba-lo_ba))

        #Overplot the data
        for dat in self.dats:
            if dat.outputs[1] is None or dat.outputs[3] is None or dat.outputs[2] is None or dat.outputs[4] is None:
                continue
            else:
                self.__datPlot(dat.outputs[0], 
                               np.mean([dat.outputs[1], dat.outputs[3]], axis=0), 
                               np.sqrt(np.power(dat.outputs[2],2)+np.power(dat.outputs[4],2)), 
                               label=dat.fhandle)
        #Overplot the bands
        if self.opBands:
            for i in range(len(self.bandLo)):
                self.__bandPlot(self.bandLo[i], self.bandHi[i])
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Transmission')
        self.__legend()
        plt.savefig('%sTransmission%s.png' % (self.saveLoc, self.fhandle))
        self.fignum += 1

    #Method for plotting output data from transmission calculation
    def plotRefl(self):
        plt.figure(self.fignum, figsize=self.figsize)
        for sim in self.sims:
            self.__simPlot(sim.outputs[0], 
                           np.mean([sim.outputs[7], sim.outputs[10]], axis=0), 
                           np.mean([sim.outputs[8], sim.outputs[11]], axis=0), 
                           np.mean([sim.outputs[9], sim.outputs[12]], axis=0), 
                           label=sim.fhandle)
            #Calculate the band-averaged reflction
            for i in range(len(self.bandLo)):
                mask = (sim.outputs[0] > self.bandLo[i])*(sim.outputs[0] < self.bandHi[i])
                mean_ba = np.trapz(np.mean([sim.outputs[7][mask], sim.outputs[10][mask]], axis=0), sim.outputs[0][mask])/(sim.outputs[0][mask][-1] - sim.outputs[0][mask][0])
                lo_ba  = np.trapz(np.mean([sim.outputs[8][mask], sim.outputs[11][mask]], axis=0), sim.outputs[0][mask])/(sim.outputs[0][mask][-1] - sim.outputs[0][mask][0])
                hi_ba  = np.trapz(np.mean([sim.outputs[9][mask], sim.outputs[12][mask]], axis=0), sim.outputs[0][mask])/(sim.outputs[0][mask][-1] - sim.outputs[0][mask][0])
                print ('%.1f GHz Band Reflection for %s = %.3f + %.3f / - %.3f' % (self.bandCenters[i], sim.fhandle, mean_ba, hi_ba-mean_ba, mean_ba-lo_ba))

        #Overplot the data
        for dat in self.dats:
            if dat.outputs[5] is None or dat.outputs[7] is None or dat.outputs[6] is None or dat.outputs[8] is None:
                continue
            else:
                self.__datPlot(dat.outputs[0], 
                               np.mean([dat.outputs[5], dat.outputs[7]], axis=0), 
                               np.sqrt(np.power(dat.outputs[6],2)+np.power(dat.outputs[8],2)), 
                               label=dat.fhandle)
        #Overplot the bands
        if self.opBands:
            for i in range(len(self.bandLo)):
                self.__bandPlot(self.bandLo[i], self.bandHi[i])
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Reflection')
        self.__legend()
        plt.savefig('%sReflection%s.png' % (self.saveLoc, self.fhandle))
        self.fignum += 1

#    def plotRefl(self, outputs, labels):
#        if not len(outputs) == len(labels):
#            raise Exception("MICROWAVE TRANSMISSION EXCEPTION: In plot.plotRefl(), number of outputs and number of labels have different lengths: %d, %d" % (len(outputs), len(labels)))
#        plt.figure(self.fignum, figsize=self.figsize)
#        for i in range(len(outputs)):
#            output = outputs[i]
#            label  = labels[i]
#            #Plot central values
#            if label.upper() in self.measHandles:
#                self.__datPlot(output[0], np.mean([output[5], output[7]], axis=0), np.sqrt(np.power(output[6],2)+np.power(output[8],2)), label=label)
#            else:
#                self.__simPlot(output[0], 
#                               np.mean([output[7], output[10]], axis=0), 
#                               np.mean([output[8], output[11]], axis=0), 
#                               np.mean([output[9], output[12]], axis=0), 
#                               #np.sqrt(np.power(output[6],2)+np.power(output[8],2)), 
#                               label=label)
#            if self.opBands:
#                if label.upper() in self.measHandles:
#                    for i in range(len(self.bandLo)):
#                        self.__bandPlot(self.bandLo[i], self.bandHi[i])
#                        mask = (output[0] > self.bandLo[i])*(output[0] < self.bandHi[i])
#                        mean_ba = np.trapz(np.mean([output[5][mask], output[7][mask]], axis=0), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        std_ba  = np.trapz(np.sqrt(np.power(output[6][mask],2)+np.power(output[8][mask],2)), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        print ('%.1f GHz Band Reflection for %s = %.3f +/- %.3f' % (self.bandCenters[i], label, mean_ba, std_ba))
#                    print ('')
#                else:
#                    for i in range(len(self.bandLo)):
#                        self.__bandPlot(self.bandLo[i], self.bandHi[i])
#                        mask = (output[0] > self.bandLo[i])*(output[0] < self.bandHi[i])
#                        mean_ba = np.trapz(np.mean([output[7][mask], output[10][mask]], axis=0), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        lo_ba  = np.trapz(np.mean([output[8][mask], output[11][mask]], axis=0), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        hi_ba  = np.trapz(np.mean([output[9][mask], output[12][mask]], axis=0), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        print ('%.1f GHz Band Transmission for %s = %.3f + %.3f / - %.3f' % (self.bandCenters[i], label, mean_ba, hi_ba-mean_ba, mean_ba-lo_ba))
#                    print ('')                    
#        plt.xlabel('Frequency [GHz]')
#        plt.ylabel('Reflection')
#        #plt.ylim([0,1])
#        self.__legend()
#        plt.savefig('%sReflection%s.png' % (self.saveLoc, self.fname))
#        self.fignum += 1

#    def plotAbso(self, outputs, labels):
#        if not len(outputs) == len(labels):
#            raise Exception("MICROWAVE TRANSMISSION EXCEPTION: In plot.plotAbso(), number of outputs and number of labels have different lengths: %d, %d" % (len(outputs), len(labels)))
#        plt.figure(self.fignum, figsize=self.figsize)
#        for i in range(len(outputs)):
#            output = outputs[i]
#            label  = labels[i]
#            #Plot central values
#            if label.upper() in self.measHandles:
#                self.__datPlot(output[0], np.mean([output[9], output[11]], axis=0), np.sqrt(np.power(output[10],2)+np.power(output[12],2)), label=label)
#            else:
#                self.__simPlot(output[0], 
#                               np.mean([output[13], output[16]], axis=0), 
#                               np.mean([output[14], output[17]], axis=0), 
#                               np.mean([output[15], output[18]], axis=0), 
#                               #np.sqrt(np.power(output[10],2)+np.power(output[12],2)), 
#                               label=label)
#            if self.opBands:
#                if label.upper() in self.measHandles:
#                    for i in range(len(self.bandLo)):
#                        self.__bandPlot(self.bandLo[i], self.bandHi[i])
#                        mask = (output[0] > self.bandLo[i])*(output[0] < self.bandHi[i])
#                        mean_ba = np.trapz(np.mean([output[9][mask], output[11][mask]], axis=0), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        std_ba  = np.trapz(np.sqrt(np.power(output[10][mask],2)+np.power(output[12][mask],2)), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        print ('%.1f GHz Band Absorption for %s = %.3f +/- %.3f' % (self.bandCenters[i], label, mean_ba, std_ba))
#                    print ('')
#                else:
#                    for i in range(len(self.bandLo)):
#                        self.__bandPlot(self.bandLo[i], self.bandHi[i])
#                        mask = (output[0] > self.bandLo[i])*(output[0] < self.bandHi[i])
#                        mean_ba = np.trapz(np.mean([output[13][mask], output[16][mask]], axis=0), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        lo_ba  = np.trapz(np.mean([output[14][mask], output[17][mask]], axis=0), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        hi_ba  = np.trapz(np.mean([output[15][mask], output[18][mask]], axis=0), output[0][mask])/(output[0][mask][-1] - output[0][mask][0])
#                        print ('%.1f GHz Band Transmission for %s = %.3f + %.3f / - %.3f' % (self.bandCenters[i], label, mean_ba, hi_ba-mean_ba, mean_ba-lo_ba))
#                    print ('')                    
#        plt.xlabel('Frequency [GHz]')
#        plt.ylabel('Absorption')
#        plt.ylim([0,1])
#        self.__legend()
#        plt.savefig('%sAbsorption%s.png' % (self.saveLoc, self.fname))
#        self.fignum += 1
        
    def plotAll(self, outputs, labels):
        self.plotTrans(output, labels)
        self.plotRefl(output,  labels)
        self.plotAbso(output,  labels)
        return

    #Generic funciton for plotting
    def __simPlot(self, freq, cent, lo, hi, label):
        p = plt.plot(freq, cent, linewidth=self.lw, label=label)
        #if not np.all(lo) == np.all(cent) and not np.all(hi) == np.all(cent):
        if not np.all(lo == 0) and not np.all(hi == 0):
            plt.fill_between(freq, lo, hi, color=p[0].get_color(), linewidth=0, alpha=0.25)
        return 1

    #Function for overplotting data from various apparatuses
    def __datPlot(self, freq, cent, err, label):
        plt.errorbar(freq, cent, yerr=err, marker='o', markersize=self.ms, linewidth=self.lw, label=label)
        return 1

    def __bandPlot(self, bandLo, bandHi):
        plt.axvspan(bandLo, bandHi, linewidth=0, color='k', alpha=0.2)
        return 1

    #Generic function for generting legend
    def __legend(self):
        plt.legend(loc='best', fontsize=self.lfsz)
        return 1

    #Function for 
