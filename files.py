

# Functions pertaining to creating, reading files and directories

# Functions contained:
#   read_control_file
#   setup_galprops
#   get_transition_info
#   setup_idcells
#   setup_mockspec
#   setup_ion_dir

import sys
import numpy as np
import os
import subprocess as sp


class runProps(object):

    '''
    A class to describe the parameters of the run of mockspec. 
    Basically holds everything read in from mockspec.config
    in read_control_file except ion information
    '''
    
    def __init__ (self):
        
        # General properties
        self.galID = 'galID'
        self.expn = '1.000'
        self.nlos = 1000
        self.maximpact = 1.5
        self.incline = 0
        self.ewcut = 0.0
        self.snr = 30
        self.ncores = 1
        self.runLoc = ''
        self.sigcellsCut = 5.
        
        # Flags
        self.runRates = 1
        self.runGenLOS = 1
        self.runCellfinder = 1
        self.runIdcells = 1
        self.runLos7 = 1
        self.runSpecsynth = 1
        self.runSysanal = 1
        self.runCullabs = 1
        self.runLocateCells = 1
        self.runSummaries = 1
        self.runTPCF = 1
        self.runPlotting = 1
        
        
class ionProps(object):

    '''
    Class to describe an ion being studied. Atributtes include the 
    ion name, the metallicity enhancement, and the instrument
    '''
    
    def __init__ (self):
        
        self.name = 'HI'
        self.xh = 0.
        self.instrument = 'COSNUV'
    
class tpcfProps(object):
    
    '''
    Class to decribe settings for TPCF
    '''

    def __init__ (self):

        self.ewLo = 0.
        self.ewHi = 10.
        self.dLo = 0.
        self.dHi = 200.
        self.azLo = 0.
        self.azHi = 90.
        self.binSize = 10.
        self.bootNum = 1000


def read_control_file():

    """
    Read in the control file, named mockspec.config
    Returns a runProps object
    """

    run = runProps() 
    tpcfP = tpcfProps()

    filename = 'mockspec.config'
    try:
        f = open(filename)
    except IOError:
        print('Missing control file: {0:s}'.format(filename))
        print('Quitting...')
        sys.exit()

    # Read past header
    f.readline()
    f.readline()

    run.galID = f.readline().split()[0]
    run.expn = f.readline().split()[0]
    run.nlos = int(f.readline().split()[0])
    run.maximpact = float(f.readline().split()[0])
    run.incline = float(f.readline().split()[0])
    run.ewcut = float(f.readline().split()[0])
    run.snr = float(f.readline().split()[0])
    run.ncores = int(f.readline().split()[0])
    run.runLoc = f.readline().split()[0]
    run.sigcellsCut = float(f.readline().split()[0])

    # Now at flags for running the various subfunctions
    f.readline()
    run.runRates = int(f.readline().split()[0])
    run.runGenLOS = int(f.readline().split()[0])
    run.runCellfinder = int(f.readline().split()[0])
    run.runIdcells = int(f.readline().split()[0])
    run.runLos7 = int(f.readline().split()[0])
    run.runSpecsynth = int(f.readline().split()[0])
    run.runSysanal = int(f.readline().split()[0])
    run.runCullabs = int(f.readline().split()[0])
    run.runLocateCells = int(f.readline().split()[0])
    run.runSummaries = int(f.readline().split()[0])
    run.runTPCF = int(f.readline().split()[0])
    run.runPlotting = int(f.readline().split()[0])

    # Now at TPCF settings
    f.readline()
    tpcfP.ewLo = float(f.readline().split()[0])
    tpcfP.hiLo = float(f.readline().split()[0])
    tpcfP.dLo = float(f.readline().split()[0])
    tpcfP.dHi = float(f.readline().split()[0])
    tpcfP.azLo = float(f.readline().split()[0])
    tpcfP.azHi = float(f.readline().split()[0])
    tpcfP.binSize = float(f.readline().split()[0])
    tpcfP.bootNum = int(f.readline().split()[0])

    # Now at ion section
    # Loop over rest of file
    ions = []
    f.readline()
    f.readline()
    for line in f:
        if line.startswith('#')==False:
            l = line.split()
            ion = ionProps()
            ion.name = l[0]
            ion.xh = l[1]
            ion.instrument = l[2]
            ions.append(ion)

    f.close()

    run.runLoc = os.getcwd()
    #props = (galID, expn, nlos, maximpact, incline, ewcut, snr, ncores, rootLoc, sigcellsCut)
    #flags = (runRates, runGenLOS, runCellfinder, runIdcells, runLos7, 
    #         runSpecsynth, runSysanal, runCullabs, runLocateCells, runSummaries, runPlotting)
    #return props, flags, ions, xh, instruments
    return run,ions,tpcfP


def setup_galprops(run, requiredLoc, summaryLoc):

    from subprocess import call

    #if not os.path.exists('gal_props.dat'):
    # Too many errors occurred with running with old gal_props
    # Now overwrite it everytime
    #command = 'cp ' +requiredLoc+ '/gal_props.dat .'
    command = 'cp {0:s}/gal_props.dat .'.format(requiredLoc)
    call(command, shell=True)

    call('cp gal_props.dat gal_props.dat.tmp', shell=True)
    f = open('gal_props.dat.tmp')
    gpf = open('gal_props.dat', 'w')

    for i in range(0,4):
        line = f.readline()
        gpf.write(line)

    gasfile = '{0:s}_GZa{1:s}.txt'.format(run.galID,run.expn)
    line = f.readline()
    gpf.write(line.replace('MW9_GZ932.a1.001.HI.txt', gasfile))
    
    line = f.readline()
    gpf.write(line.replace('MW9', run.galID))

    line = f.readline()
    gpf.write(line.replace('1.001', run.expn))

    gpf.write(f.readline())
    gpf.write(f.readline())
    
    gpf.write(summaryLoc)
    f.close()
    gpf.close()

    call('rm gal_props.dat.tmp', shell=True)



def get_transition_info(ion, codeLoc):

    # Returns the atomic number and existion level 
    # of this ion
    # Open the Mockspec.transistions file
    f = open(codeLoc+'/data/Mockspec.transitions')

    for line in f:
        fion = line.split()[4]
        if fion==ion:
            element = line.split()[1]
            k = line.split()[2]
            j = line.split()[3]
            f.close()
            return element, k, j




def setup_idcells(gasfile, ion_list):
    
    from subprocess import call

    if not os.path.exists('idcells_props.dat'):
        command = 'cp /home/matrix3/jrvander/requiredfiles/idcells_props.dat .'
        call(command, shell=True)

    call('cp idcells_props.dat idcells_props.tmp', shell=True)
    tmpf = open('idcells_props.tmp')
    datf = open('idcells_props.dat', 'w')

    for i in range(0,3):
        line = tmpf.readline()
        datf.write(line)

    line = tmpf.readline()
    datf.write(line.replace('dwarf9o_GZa1.001.txt', gasfile))
    
    line = tmpf.readline()
    datf.write(line)

    line = tmpf.readline()
    ion_string = ''
    for ion in ion_list:
        ion_string = ion_string+ion+' '
    datf.write(ion_string + '\t\t# Ions to use')
    tmpf.close()
    datf.close()

    call('rm idcells_props.tmp', shell=True)






#def setup_mockspec(ion_list, instr_list, ewcut, snr, xh_list, requiredLoc):
def setup_mockspec(ions,run,requiredLoc):

    '''
    Ensures that the appropriate Mockspc control files are present in 
    the inclination directory.  Starts in the inclination directory. 
    Copies Mockspec.runpars, Mockspec.instruments, and Mockspec.transitions 
    from the control directory if they are not already present.
    '''

    from subprocess import call

    cwd = os.getcwd()
    runfile = '{0:s}/Mockspec.runpars'.format(cwd)
    instfile = '{0:s}/Mockspec.instruments'.format(cwd)
    transfile = '{0:s}/Mockspec.transitions'.format(cwd)

    # Need Mockspec.runpars, Mockspec.instruments, Mockspec.transitions
    basecommand = 'cp '+requiredLoc
    if not os.path.isfile(runfile):
        command = basecommand + 'Mockspec.runpars .'
        call(command, shell=True)

    if not os.path.isfile(instfile):
        command = basecommand + 'Mockspec.instruments .'
        call(command, shell=True)

    if not os.path.isfile(transfile):
        command = basecommand + 'Mockspec.transitions .'
        call(command, shell=True)


    # Alter Mockspec.runpars
    # The others should never be altered
    call('cp Mockspec.runpars Mockspec.tmp', shell=True)
    tmpf = open('Mockspec.tmp')
    datf = open('Mockspec.runpars', 'w')

    line = tmpf.readline()
    datf.write(line)

    slev = '5.0'  # Significance level of detection
    
    for ion in ions:

        element, Z, stage = get_transition_info(ion.name, requiredLoc)
        if element =='hydrogen':
            vmax = '10000.'
        else:
            vmax = '1000.'
        element = '\''+element+'\''
        line = ('{0:<11s} {1:<5s} {2:<6s} {3:<7s} '
                '{4:<6f} {5:<5f} {6:<6s} {7:<s}\n'.format(element, 
                stage, ion.xh, ion.instrument, run.ewcut, 
                run.snr, vmax, slev))

        datf.write(line)
        
    call('rm Mockspec.tmp', shell=True)
    tmpf.close()
    datf.close()


def get_transition_info(ion, requiredLoc):

    """
    Returns the atomic number and exiciation level
    of this ion as defined in Mockspec.transitions
    """
    filename = requiredLoc+'/Mockspec.transitions'
    f = open(filename)
    for line in f:
        fion = line.split()[4]
        if fion==ion:
            element = line.split()[1]
            k = line.split()[2]
            j = line.split()[3]
            f.close()
            return element, k, j
    


#def setup_ion_dir(ion, galID, expn, codeLoc):
def setup_ion_dir(ion, run, codeLoc):
    
    """ 
    Creates a directory for this ion if one does not exit
    Copies the all required files into this directory
    Returns the path to the ion directory
    """
    ionloc = '{0:s}/i{1:d}/{2:s}'.format(run.runLoc,int(run.incline),ion.name)
    
    if not os.path.exists(ionloc):
        command = 'mkdir '+ion.name
        try:
            sp.check_call(command, shell=True)
        except:
            print('Could not complete {0:s}'.format(command))

    # Create the los.list file
    command = ('ls *{0:s}.los*dat > qso.list && mv qso.list '
                './{0:s}/'.format(ion.name))
    try:
        sp.check_call(command, shell=True, stderr=sp.PIPE)
    except:
        print('Could not complete {0:s}'.format(command))
    
    # Copy the cell files, ion boxes and the lines files into the ion directory
    ionbox = '{0:s}_GZa{1:s}.{2:s}.h5'.format(run.galID,run.expn,ion.name)
    #command = 'cp '+ionbox+' ./'+ion+'/'
    command = 'cp {0:s} ./{1:s}'.format(ionbox,ion.name)
    try:
        sp.call(command, shell=True)
    except: 
        print('Could not complete {0:s}'.format(command))

    #command = 'mv *'+ion+'.los*.dat ./'+ion+'/'
    command = 'mv *{0:s}.los*.dat ./{0:s}/'.format(ion.name)
    try:
        sp.check_call(command, shell=True)
    except:
        print('Could not complete {0:s}'.format(command))
        
    command = 'cp lines* ./'+ion.name+'/'
    try:
        sp.check_call(command, shell=True)
    except:
        print('Could not complete {0:s}'.format(command))
   

    # Copy the Mockspec files from the parent directory to here
    # New way: pull files from parent directory, not controls
    command = 'cp ./Mockspec.instruments ./{1:s}/'.format(codeLoc,ion.name)
    try:
        sp.check_call(command, shell=True)
    except:
        print('Could not complete {0:s}'.format(command))

    command = 'cp ./Mockspec.transitions ./{1:s}/'.format(codeLoc,ion.name)
    try:
        sp.check_call(command, shell=True)
        print('Copying Mockspec.transitions: \n{0:s}'.format(command))
    except:
        print('Could not complete {0:s}'.format(command))

    command = 'cp ./Mockspec.runpars ./{1:s}/'.format(codeLoc,ion.name)
    try:
        sp.check_call(command, shell=True)
    except:
        print('Could not complete {0:s}'.format(command))

    # Old way: Copy from controls
    # This overwrites any local changes made, such as to the transitions file to
    # observe specific transitions
#    command = 'cp {0:s}/controls/Mockspec.instruments ./{1:s}/'.format(codeLoc,ion)
#    try:
#        sp.check_call(command, shell=True)
#    except:
#        print 'Could not complete {0:s}'.format(command)
#
#    command = 'cp {0:s}/controls/Mockspec.transitions ./{1:s}/'.format(codeLoc,ion)
#    try:
#        sp.check_call(command, shell=True)
#    except:
#        print 'Could not complete {0:s}'.format(command)
#
#    command = 'cp ./Mockspec.runpars ./{1:s}/'.format(codeLoc,ion)
#    try:
#        sp.check_call(command, shell=True)
#    except:
#        print 'Could not complete {0:s}'.format(command)

    # Check that nothing exists in the parent directory
    #command = 'ls *.{0:s}.*los*dat | wc -l'.format(ion.name)
    command = 'ls *{0:s}*los*dat &> dum.txt && grep -c dat$ dum.txt'.format(ion.name)
    try:
        numDatFiles = int(sp.check_output(command, shell=True).strip())
    except sp.CalledProcessError:
        numDatFiles = 0
        pass

    if numDatFiles!=0:
        print('ERROR in setup_ion_dir in files.py')
        print('Problem moving .dat file for {0:s}'.format(ion.name))
        print('Exitting....')
        sys.exit()
    return ionloc


#def setup_inclination_dir(incline, ions, runRates, galID, expn):
def setup_inclination_dir(run, ions):

    """
    Sets up a directort called i<incline> in the expansion factor
    directory. Moves all rates output files into it if rates
    was run. Returns the file path to the new directory
    """
    inc = int(run.incline)    

    # Check to see if the inclination directory already
    # exits
    if not os.path.exists('./i{0:d}'.format(inc)):
        command = 'mkdir i{0:d}'.format(inc)
        sp.call(command, shell=True)

    cwd = os.getcwd()
    incLoc = '{0:s}/i{1:d}/'.format(cwd, inc)

    # Check to see if the ion boxes are already in the directory
    for ion in ions:
        boxName = '{0:s}_GZa{1:s}.{2:s}.h5'.format(run.galID, run.expn, ion.name)
        if not os.path.isfile('{0:s}/{1:s}'.format(incLoc,boxName)):
            command = 'cp {0:s} {1:s}/'.format(boxName, incLoc)
            try:
                sp.check_call(command, shell=True)
            except sp.CalledProcessError:
                print('Error in files.py in setup_inclination_dir')
                print('Could not copy {0:s} into {1:s}'.format(boxName, incLoc))
                continue
                
    # Copy the rest of the control files into the directory
    boxname = '{0:s}_GZa{1:s}.txt'.format(run.galID, run.expn)
    filenames = ['mockspec.config', 'gal_props.dat', 'galaxy.props', boxname]
    for fn in filenames:
        command = 'cp {0:s} ./i{1:d}/'.format(fn, inc)
        try:
            sp.check_call(command, shell=True)
        except sp.CalledProcessError:
            print('Cannot find {0:s} in setup_inclination_dir'.format(fn))
            sys.exit()


    # Return the path to the new directory
    return incLoc


#def setup_galaxy_props(summaryLoc, galID, expn, inc):
def setup_galaxy_props(run, sumFile):
        
    """
    Creates a file called galaxy.props in the root directory.
    Contains relevent galaxy information from the <galID>.dat file
    stored in the codeLoc's summary directory
    """
    
#    cwd = os.getcwd()
    # Read in properties from the summary file
#    summaryFile = '{0:s}.dat'.format(galID)
#    try:
#        f = open(summaryLoc+summaryFile)
#        summ = 1
#    except IOError:
#        try:
#            f = open('../output/rotmat_a{0:s}.txt'.format(expn))
#            summ = 2
#        except IOError:
#            print 'Could not open summary file:\n\t{0:s}{1:s}'.format(summaryLoc, summarFile)
#            print 'Coult not find rotmat file'
#            sys.exit()

    try:
        f = open(sumFile, 'r')
    except IOError:
        print('Cannot open {0:s} in setup_galaxy_props'.format(sumFile))
        print('Exitting...')
        sys.exit()
    f.readline()
    l = f.readline().split()
    a = '{0:.3f}'.format(float(l[0]))
    if a==run.expn:
        found = 1
        redshift = float(l[1])
        mvir = float(l[2])
        rvir = float(l[3])
    else:
        print('Could not find {0:s} in {1:s}'.format(run.expn, sumFile))
        sys.exit()
    f.close()

    # Write the galaxy.props file
    with open('galaxy.props', 'w') as f:

        s = 'galID          {0:s}\n'.format(run.galID)
        f.write(s)
        s = 'Expn           {0:s}\n'.format(run.expn)
        f.write(s)
        s = 'Redshift       {0:.3f}\n'.format(redshift)
        f.write(s)
        s = 'Mvir           {0:.4e}\n'.format(mvir)
        f.write(s)
        s = 'Rvir           {0:.4f}\n'.format(rvir)
        f.write(s)
        s = 'Inclination    {0:.1f}\n'.format(run.incline)
        f.write(s)



    


#def rename(galID, expn, ion, incline, runLocateCells, runCullabs):
def rename(run,ion):

    '''
    Renames files output by the pipeline to incline the 
    inclination angle
    '''
    inc = int(run.incline)
    rootLoc = '{0:s}/i{1:d}/{2:s}/'.format(run.runLoc,
                int(run.incline),ion.name)
    if run.runCullabs==1:
        # Need to rename ALL.sysabs file
        allName = '{0:s}.{1:s}.a{2:s}.ALL.sysabs.h5'.format(run.galID,
                     ion.name,run.expn)
        newName = allName.replace('ALL','i{0:d}.ALL'.format(int(run.incline)))
        command = 'mv {0:s}{1:s} {0:s}{2:s}'.format(rootLoc,
                    allName, newName)
        try:
            sp.check_call(command,shell=True,stderr=sp.PIPE)
        except:
            print('Error in rename running \n\t{0:s}'.format(command))
            print('\tCWD: ',os.getcwd())
       

def timing_setup(startTime,run,ions,tpcfProp):

    '''
    Generates a file containing the timing of the run
    '''

    import os
    import time
    import itertools as it
    import datetime as dt
    import textwrap as tw

    fname = 'timing.out'
    f = open(fname,'w')

    # Print starting time
    sTime = time.strftime('%H:%M:%S',time.localtime(startTime))
    sDate = time.strftime('%d-%m-%Y',time.localtime(startTime))
    line = '\nStart Time = {0:s}\tStart Date = {1:s}\n\n'.format(sTime,sDate)
    f.write(line)

    # Print run params
    line = 'Run Settings:\n'
    line += 'galID = {0:<10s}\ta = {1:<10s}\n'.format(run.galID,run.expn)
    line += 'nLOS  = {0:<10d}\ti = {1:<10.1f}\n'.format(run.nlos,run.incline)
    line += 'ewCut = {0:<10.1f}\tSNR = {1:<10.1f}\n'.format(run.ewcut,run.snr)
    line += 'ncores = {0:<10d}\tsigCut = {1:<10.1f}\n\n'.format(run.ncores,run.sigcellsCut)
    f.write(line)

    # Print out flags
    flags = [run.runRates,run.runGenLOS,run.runCellfinder,
             run.runIdcells,run.runLos7,run.runSpecsynth,
             run.runSysanal,run.runCullabs,run.runLocateCells,
             run.runSummaries,run.runTPCF,run.runPlotting]
    flags = [True if i==1 else False for i in flags]
    functions = ['rates','genLOS','cellfinder','idCells','los7',
                 'specsynth','sysanal','cullabs','locatecells',
                 'summaries','tpcf','plotting']
    runFuncs = list(it.compress(functions,flags))

    line = 'Functions Run:\n'
    line += tw.fill(', '.join(runFuncs))
    f.write(line+'\n\n')


    # Print out TPCF settings if TPCF was performed
    if run.runTPCF==1:
        line = 'TPCF Settings:\n'
        line += 'EW limits = {0:.1f} - {1:.2f}\n'.format(tpcfProp.ewLo,tpcfProp.ewHi)
        line += 'Impact limits = {0:.1f} - {1:.2f}\n'.format(tpcfProp.dLo,tpcfProp.dHi)
        line += 'Az limits = {0:.1f} - {1:.2f}\n'.format(tpcfProp.azLo,tpcfProp.azHi)
        line += 'Binsize = {0:.1f} \tnBoot = {1:d}\n\n'.format(tpcfProp.binSize,tpcfProp.bootNum)
        f.write(line)

    # Print ions:
    line = 'Ions:\n'
    ionLine = ['{0:s} ({1:s})'.format(i.name,i.instrument) for i in ions]
    line += ', '.join(ionLine)
    f.write(line)
   
    # Start timing information
    timestring = '{0[0]:<14s}\t{0[1]:<10s}\t{0[2]:<8s}\t{0[3]:<12s}\n'
    header = 'Event Date Time Elapsed'.split()
    f.write('\n\n')
    f.write(timestring.format(header))
    
    # Print out starting times
    current = time.time()
    dateString = time.strftime('%d-%m-%Y',time.localtime(current))
    timeString = time.strftime('%H:%M:%S',time.localtime(current))
    elapsed = str(dt.timedelta(seconds=current-startTime))
    line = ['Starting',dateString,timeString,elapsed]
    f.write(timestring.format(line))

    fullFileName = os.path.realpath(f.name)
    f.close()
    
    print(fullFileName,flush=True)
    return fullFileName,timestring



def record_time(timef,timestring,event,current,startTime):
    '''
    Records the current time in to the timing file
    '''

    import time
    import datetime as dt

    with open(timef,'a') as f:
        dateString = time.strftime('%d-%m-%Y',time.localtime(current))
        timeString = time.strftime('%H:%M:%S',time.localtime(current))
        elapsed = str(dt.timedelta(seconds=current-startTime))
        line = [event,dateString,timeString,elapsed]
        f.write(timestring.format(line))
    
    












