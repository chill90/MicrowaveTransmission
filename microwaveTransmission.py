#Using python 2.7.2
import sys   as sy
import          os
import numpy as np

#Custom classes
import src.simulate    as sm
import src.plot        as pl
import src.measurement as ms

#Default configuration files
layerFileDef = ('config'+os.sep+'layers'+os.sep+'layers.txt')
simFileDef   = ('config'+os.sep+'simulation'+os.sep+'simInputs.txt')

#Default save location
saveLoc = os.path.abspath('Data')
plotLoc = os.path.abspath('Plots')

#Default save header
saveHdr = "%-11s%-26s%-26s%-26s%-26s%-26s%-26s" % ("Freq [GHz]", "P Trans (mean +/- std)", "S Trans (mean +/- std)", "P Refl (mean +/- std)", "S Refl (mean +/- std)", "P Absorb (mean +/- std)", "S Absorb (mean +/- std)")

#For now, only ability is to simulate using Hou code
allowedCmds = ['LF', 'SF', 'DT']
allowedInst = ['UMich_Reflectometer',
               'Dick_CoherentSource']
def help(val):
    print ("\nERROR: could not understand '%s'" % (val))
    print ("Usage: python microwaveTransmission.py -lf [layerFile] -sf [simFile] -nt [1] -sm [HOU]")
    print ("-lf: layer file(s) that contains the dielectric layer parameters. Default value = %s" % (layerFileDef))
    print ("-dt: data file(s) that contain tranmission/reflection/absorption vs frequency. File must contain indicator of setup = %s" % (','.join(allowedInst)))
    print ("-sf: file that contains the simulation inputs. Default value = %s" % (simFileDef))
    print ("Allowed simulation methods are:")
    print ("'HOU': uses matrix formalism laid out in Hou et al. Not suitable for birefringent stacks")
    sy.exit()

#Collect command-line arguments
argstr = ' '.join(sy.argv[1:])
args = argstr.split('-')[1:]
layerFiles = []
dataFiles  = []
simFile    = simFileDef
for arg in args:
    cmd = arg.split()[0]; vals = list(arg.split()[1:])
    if cmd.upper() not in allowedCmds:
        help(cmd)
    if cmd.upper() == 'LF':
        layerFiles = vals
        for layerFile in layerFiles:
            if not os.path.isfile(layerFile):
                sy.exit("\nERROR: could not find layer file '%s'\n" % (layerFile))
    elif cmd.upper() == 'SM':
        simFile = val
        if not os.path.isfile(simFile):
            sy.exit("\nERROR: could not find sim file '%s'\n" % (simFile))
    elif cmd.upper() == 'DT':
        dataFiles = vals
        for dataFile in dataFiles:
            if not os.path.isfile(dataFile):
                sy.exit("\nERROR: could not find data file '%s'\n" % (dataFiles))
            inst = None
            for instrument in allowedInst:
                if instrument in dataFile:
                    inst = instrument
            if inst is None:
                print ("\nERROR: could not indenfity allowed instrument in passed file '%s' for overplotting" % (dataFile))
                help(cmd)

#Generate and execute simulation and plotting objects
sims = [sm.Simulate(layerFile=os.path.abspath(layerFile), simFile=os.path.abspath(simFile)) for layerFile in layerFiles]
for sim in sims:
    sim.calc()

#Write simulated output to a text file
fhandles = []
for i in range(len(layerFiles)):
    layerFile = layerFiles[i]
    output    = sims[i].outputs
    fhandle = layerFile.split('.')[0].split('/')[-1]
    fhandles.append(fhandle)
    fname = ('%s%ssimOutput_%s.txt' % (saveLoc, os.sep, fhandle))
    np.savetxt(fname, np.array(output).T, fmt="%-12.4f", header=saveHdr)

#Gather measured data
dats = [ms.Measurement(dataFile=os.path.abspath(dataFile)) for dataFile in dataFiles]
for dat in dats:
    dat.loadData()
#fhandles = fhandles + [dataFile.split('.')[0].split('/')[-1] for dataFile in dataFiles]

#Plot all data
plt = pl.Plot(sims, dats, saveLoc=plotLoc)
plt.plotTrans()
plt.plotRefl()
