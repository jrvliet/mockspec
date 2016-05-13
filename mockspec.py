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
import ratesControl as rc
import funcs.locatecells.locatecells as lc
import funcs.idcells.idcells as ic
import runFunctions as rf 
import funcs.sigcells.sigcells as sc
import funcs.analysis.analysis_control as ac

# Get the location where the code lives
pathname = os.path.dirname(sys.argv[0])
codeLoc = os.path.abspath(pathname)
summaryLoc = codeLoc+'/summaries/'
requiredLoc = codeLoc+'/controls/'

#  Read in the control file
print '\n\nReading in control file...'
props, flags, ions, xh, instruments = fi.read_control_file()
galID, expn, nlos, maximpact, incline, ewcut, snr, ncores, rootLoc, sigcellsCut = props
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

# Test the summary location
sumFile = '{0:s}/rotmat_a{1:s}.txt'.format(os.getcwd(), expn)

#sumFile = '{0:s}/{1:s}.dat'.format(summaryLoc, galID)
#if not os.path.isfile(sumFile):
    # Summary file does not exit
    # make it
#    rotfile = '../output/rotmat_a{0:s}.txt'.format(expn)
##    f = open(sumFile,'w')
#    frot = open(rotfile)
#    header = frot.readline()
#    print header
#    f.write(header)
#    f.write(header)
#    for line in frot:
#        print line
#        f.write(line)
#    frot.close()
#    f.close()



# Generate gal_props.dat file
print '\n\nGenerating gal_props.dat...'
fi.setup_galprops( galID, expn, requiredLoc, summaryLoc )

# Generate galaxy.props file, needed for analysis codes
print '\n\nGenerating galaxy.props...'
fi.setup_galaxy_props(sumFile, galID, expn, incline)

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
#  Set up inclination directory
#
#####
mainLoc = os.getcwd()
incLoc = fi.setup_inclination_dir(incline, ions, runRates, galID, expn)
os.chdir(incLoc)



##### 
#  
#  Generate the lines of sight
#
#####
if runGenLOS==1:
    print '\nGenerating lines of sight...'
    rf.genLOS( codeLoc, galID, mainLoc, expn, incline, 
               nlos, maximpact, ncores)   
else:
    print 'Skipping genLOS...'

##### 
#  
#  Run lines of sight through box with cellfinder
#
#####
if runCellfinder==1:
    print '\nRunning LOS through box...'
    rf.runCellfinder(codeLoc, ncores)
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

# Setup Mockspec 
fi.setup_mockspec(ions, instruments, ewcut, snr, xh, requiredLoc)


# Start looping over ions
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
        lc.locateSigCells(galID, expn, ion, sigcellsCut, codeLoc, incline)
#        lc.sigCells(galID, expn, ion, sigcellsCut, codeLoc)
    else:
        print '\tSkipping locatecells...'

    
    ##### 
    #  
    #  Rename files
    #  Add inclination angle to filenames
    #
    #####
    print '\n\t Renaming files...'
    if runCullabs==1 or runLocateCells==1:
        fi.rename(galID, expn, ion, incline, runLocateCells, runCullabs)
    else:
        print 'What was the point...'


    # Move back up to the parent directory
    os.chdir('..')

    
# Create basic plots
if runAnalysis==1:
    print '\n\nGenerating plots'
    ac.make_plots(ions)


