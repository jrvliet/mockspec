#
#  mockspec.py
# 
#  Full pipeline for performing synthetic QSO observations
#  of ART simulations of galaxy CGM

# General libraries
import numpy as np
import tables as tb
import os
import pandas
import subprocess as sp
import sys

# Code specific libraries
import files as fi
import genLOS as gl
import ratesControl as rc
import funcs.locatecells.locatecells as lc
import funcs.idcells.idcells as ic
import runFunctions as rf 

pathname = os.path.dirname(sys.argv[0])
codeLoc = os.path.abspath(pathname)
summaryLoc = codeLoc+'/summaries/'
requiredLoc = codeLoc+'/comtrols/'

#  Read in the control file
print '\n\nReading in control file...'
props, flags, ions, xh, instruments = fi.read_control_file()
galID, expn, nlos, maximpact, incline, ewcut, snr, ncores, rootLoc, sigcellsCut= props
runRates, runGenLOS, runCellfinder, runIdcells, runLos7, runSpecsynth, runSysanal, runCullabs, runLocateCells = flags

# Genearte the name of the gasfile
gasfile = galID+'_GZa'+expn+'.txt'
print '\nRun Parameters:'
print '\tGasfile:      ', gasfile
print '\tIons:         ', ions
print '\tNLOS:         ', nlos
print '\tMax impact:   ', maximpact
print '\tIncline:      ', incline
print '\tEWCut:        ', ewcut
print '\tSNR:          ', snr
print '\tNCores:       ', ncores
print '\tRoot Loc:     ', rootLoc
print '\tSigcells Cut: {0:.0%}'.format( sigcellsCut/100. )

print '\nRun Flags:'
print '\tRates:       {0:d}'.format(runRates)
print '\tGenLOS:      {0:d}'.format(runGenLOS)
print '\tCellfinder:  {0:d}'.format(runCellfinder)
print '\tIDCells:     {0:d}'.format(runIdcells)
print '\tLos7:        {0:d}'.format(runLos7)
print '\tSpecsynth:   {0:d}'.format(runSpecsynth)
print '\tSysanal:     {0:d}'.format(runSysanal)
print '\tCullabs:     {0:d}'.format(runCullabs)
print '\tLocateCells: {0:d}'.format(runLocateCells)




# Generate gal_props.dat file
print '\n\nGenerating gal_props.dat...'
fi.setup_galprops( galID, expn, requiredLoc, summaryLoc )


##### 
#  
#  Run rates
#
#####
if runRates==1:
    print '\nRates:'
    print '\t Setting up rates.inp...'
    rc.setup_rates_control( gasfile, expn, ions, requiredLoc)
    print '\t Setting up rates.outfiles...'
    rc.setup_rates_outputs(galID, expn, ions, codeLoc, requiredLoc) 
    print '\t Running rates...'
    rc.run_rates(codeLoc)
else:
    print 'Skipping rates'

##### 
#  
#  Generate the lines of sight
#
#####
if runGenLOS==1:
    print '\nGenerating lines of sight...'
    gl.genLines(galID, gasfile, summaryLoc, expn, incline, nlos, maximpact, ncores)
else: 
    print 'Skipping genLOS...'

##### 
#  
#  Run lines of sight through box with cellfinder
#
#####
if runCellfinder==1:
    print '\nRunning LOS through box...'
    gl.runCellfinder(codeLoc)
else: 
    print 'Skipping cellfinder...'


##### 
#  
#  Identify cells in the ion boxes
#
#####
if runIdcells==1:
    print '\nIdentifying probed cell properites...'
    ic.idcells(galID, expn, ions, codeLoc)
else:
    print 'Skipping IDcells...'



#Start looping over ions
print '\nBegin looping over ions...'
for ion in ions:

    print '\rIon: ', ion

    # Setup the ion directory
    ionloc = fi.setup_ion_dir(ion, galID, expn, codeLoc) 

    # Move into the ion directory
    os.chdir(ionloc)

    ##### 
    #  
    #  los7
    #  Determine path length through each cell and create the
    #  .losdata and .lines files
    #
    #####
    if runLos7==1:
        print '\n\tDetermining cell path length and applying rough cut'
        rf.los7(codeLoc)
    else:
        print '\tSkipping los7...'


    ##### 
    #  
    #  specsynth
    #  Generate the synthetic spectra
    #
    #####
    if runSpecsynth==1:
        print '\n\tGenerating spectra'
        rf.specsynth(codeLoc)
    else:
        print '\tSkipping specsynth...'

    ##### 
    #  
    #  sysanal
    #  Analyze the spectra
    #
    #####
    if runSysanal==1:
        print '\n\tAnalyzing spectra...'
        rf.sysanal(codeLoc)
    else:
        print '\tSkipping sysanal...'


    ##### 
    #  
    #  cullabs
    #  Create sysabs file
    #
    #####
    if runCullabs==1:
        print '\n\tCreating sysabs'
        rf.cullabs(codeLoc)
    else:
        print '\tSkipping sysabs...'

    

    ##### 
    #  
    #  locatecells
    #  Identify the cells that are significant contributers to the absorption
    #
    #####
    if runLocateCells==1:
        print '\n\tIdentifying significant cells'
        lc.locateSigCells(galID, expn, ion, sigcellsCut, codeLoc)
    else:
        print '\tSkipping locatecells...'
    os.chdir('..')

    



