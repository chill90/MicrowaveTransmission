#Using python 2.7.2
import sys   as sy
import          os
import numpy as np

#Custom classes
import src.simulate as sm
import src.plot     as pl

#Default configuration files
layerFileDef = ('config'+os.sep+'layers'+os.sep+'layers.txt')
simFileDef   = ('config'+os.sep+'simulation'+os.sep+'simInputs.txt')

#Default save location
saveLoc = os.path.abspath('Data')
plotLoc = os.path.abspath('Plots')

#Default save header
saveHdr = "%-11s%-26s%-26s%-26s%-26s%-26s%-26s" % ("Freq [GHz]", "P Trans (mean +/- std)", "S Trans (mean +/- std)", "P Refl (mean +/- std)", "S Refl (mean +/- std)", "P Absorb (mean +/- std)", "S Absorb (mean +/- std)")

#For now, only ability is to simulate using Hou code
allowedCmds = ['LF', 'SF']
allowedVals = ['HOU']
def help(val):
    print ("\nERROR: could not understand '%s'" % (val))
    print ("Usage: python microwaveTransmission.py -lf [layerFile] -sf [simFile] -nt [1] -sm [HOU]")
    print ("-lf: file that contains the dielectric layer parameters. Default value = %s" % (layerFileDef))
    print ("-sf: file that contains the simulation inputs. Default value = %s" % (simFileDef))
    print ("Allowed simulation methods are:")
    print ("'HOU': uses matrix formalism laid out in Hou et al. Not suitable for birefringent stacks")
    sy.exit()

#Collect command-line arguments
argstr = ' '.join(sy.argv[1:])
args = argstr.split('-')[1:]
layerFile = layerFileDef
simFile   = simFileDef
for arg in args:
    cmd, val = arg.split()
    if cmd.upper() not in allowedCmds:
        help(cmd)
    if cmd.upper() == 'LF':
        layerFile = str(val)
        if not os.path.isfile(layerFile):
            sy.exit("\nERROR: could not find layer file '%s'\n" % (layerFile))
    elif cmd.upper() == 'SM':
        simFile = str(val)
        if not os.path.isfile(layerFile):
            sy.exit("\nERROR: could not find sim file '%s'\n" % (layerFile))

#Calculate transmission
sim = sm.Simulate(layerFile=os.path.abspath(layerFile), simFile=os.path.abspath(simFile))
output = sim.calc()
#Write output to a text file
fhandle = layerFile.split('.')[0].split('/')[-1]
if fhandle == '':
    fname = ('%s%ssimOutput.txt' % (saveLoc, os.sep))
else:
    fname = ('%s%ssimOutput_%s.txt' % (saveLoc, os.sep, fhandle))
np.savetxt(fname, np.array(output).T, fmt="%-12.4f", header=saveHdr)

#Plot transmission
plt = pl.Plot(saveLoc=plotLoc, fhandle=fhandle)
plt.plotAll(output)
