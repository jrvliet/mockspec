
'''
Codes to convert the ALL.sysabs files and 
abs_cells files to HDF5 format for
'''

import numpy as np
import pandas as pd
import glob
import losSummary as ls

def az(phi):

    '''
    Converts the position angle to azimuthal angle,
    the angle between the LOS and the closes major axis
    '''
    if phi<90:
        return phi
    elif phi<180:
        return 180.-phi
    elif phi<270:
        return phi-180
    else:
        return 360-phi
    

def sysabs_to_hdf5(run,ion,codeLoc):

    '''
    Converts the ALL.sysabs file into the HDF5 format
    '''
    
    # Get the name of the file 
    allfile = '{0:s}/i{1:d}/{2:s}/{3:s}.{2:s}.a{4:s}.ALL.sysabs'.format(
                run.runLoc,int(run.incline),ion.name,run.galID,run.expn)
    #allfile = (glob.glob('*.ALL.sysabs'))[0]

    # Create the name of the HDF5 file
    hdf5file = allfile+'.h5'
    hdf5file = hdf5file.replace('ALL','i{0:d}.ALL'.format(int(run.incline)))

    # Get the header
    header = ['los','D','zabs','v-','v+','EW_r','dEW_r','DR','dDR',
               'SL','Vbar','dVbar','Vsprd','dVsprd','Vasym','dVasym',
                'lgt','dtau-','dtau+','logN','dNcol-','dNcol+']

    # Read in the data
    try:
        #data = np.loadtxt(allfile, skiprows=1)
        data = pd.read_csv(allfile,sep='\s+',names=header,skiprows=1)
    except (ValueError,IOError) as e:
        print 'Error in sysabs_to_hdf5 in mkHDF5.py'
        print 'while reading in {0:s}'.format(allfile)
        return 1
    
    # Read in lines.info to get the azimuthal angle
    linesfile = 'lines.info'
    phi = np.loadtxt(linesfile, skiprows=2, usecols=(2,), unpack=True)
    
    # Insert into the sysabs data file
    #data = np.insert(data, 2, phi, axis=1)
    data['phi'] = phi
    data['azimuthal'] = data['phi'].apply(az)
    #header.insert(2, 'phi')
    

    # WRite data to HDF file
    #df = pd.DataFrame(data, columns=header)
    data.to_hdf(hdf5file, 'data', mode='w')


def regabs_to_hdf5(run,ion,codeLoc):

    '''
    Converts the ALL.regabs file into the HDF5 format
    '''
    
    # Get the name of the file 
    allfile = '{0:s}/i{1:d}/{2:s}/{3:s}.{2:s}.a{4:s}.ALL.regabs'.format(
                run.runLoc,int(run.incline),ion.name,run.galID,run.expn)
    #allfile = (glob.glob('*.ALL.regabs'))[0]

    # Create the name of the HDF5 file
    hdf5file = allfile + '.h5'

    # Get the header
    header = ['los','D','reg','zabs','v-','v+','EW_r','dEW_r','DR',
                'dDr','SL','Vbar','dVbar','Vsprd','dVsprd',
                'Vasym','dVasym','dum','ytick']

    # Read in the data
    try:
        #data = np.loadtxt(allfile, skiprows=1)
        data = pd.read_csv(allfile,sep='\s+',skiprows=1,names=header)
    except (ValueError,IOError) as e:
        print 'Error in sysabs_to_hdf5 in mkHDF5.py'
        print 'while reading in {0:s}'.format(allfile)
        return 1
    
    # Read in lines.info to get the azimuthal angle
    linesfile = 'lines.info'
    angle = np.loadtxt(linesfile, skiprows=2, usecols=(2,), unpack=True)
    
    # Insert into the sysabs data file
    #data = np.insert(data, 2, phi, axis=1)
    #header.insert(2, 'phi')
    
    phi = []
    for i in range(len(data)):
        los = int(data['los'].iloc[i])
        phi.append(angle[los-1])
    data['phi'] = phi
    
    # WRite data to HDF file
    #df = pd.DataFrame(data, columns=header)
    data.to_hdf(hdf5file, 'data', mode='w')


def abscells_to_hdf5(codeLoc):

    '''
    Converts the abs_cells.dat file into the HDF5 format
    '''
    
    # Get the name of the file 
    filename = (glob.glob('*.abs_cells.dat'))[0]

    # Create the name of the HDF5 file
    hdf5file = filename.replace('dat','h5')

    print filename
    # Get the header
    with open(filename) as f:
        header = f.readline().strip().split()
    print header
    print len(header)

    # Read in the data
    data = np.loadtxt(filename, skiprows=1)
    print data.shape    

    # WRite data to HDF file
    df = pd.DataFrame(data, columns=header)
    df.to_hdf(hdf5file, 'data', mode='w')


def gasbox_to_hdf5(codeLoc, ions, run):

    '''
    Converts the ion boxes file into the HDF5 format
    '''

    # Convert the ion boxes
    for ion in ions:
        
        filename = '{0:s}_GZa{1:s}.{2:s}.txt'.format(run.galID,
                                                run.expn,ion.name)
        header = ['cell_size', 'x', 'y', 'z', 'vx', 'vy', 'vz',
                        'nH', 'temperature', 'SNII', 'SNIa', 'nAtom',
                        'fIon', 'nIon', 'alpha_sol', 'alpha_Zmet',
                        'ID', 't_ph', 't_rec', 't_coll', 't_cool']
        
        box_conversion(filename,header)
        
    # Convert the normal gas box
    filename = '{0:s}_GZa{1:s}.txt'.format(run.galID,run.expn) 
    header = ['cell_size', 'x', 'y', 'z', 'vx', 'vy', 'vz', 
                        'density', 'temperature', 'SNII', 'SNIa']
    box_conversion(filename,header)

    # Convert the timescale gas box
    filename = '{0:s}_GZa{1:s}.tcdat.txt'.format(run.galID,run.expn)
    header = ['ID', 'z', 'nH', 'T', 'Zcell', 'ne', 'ntot', 
                      'nH/nH', 'nHe/nH', 'nC/nH', 'nN/nH', 'nO/nH', 
                      'nNe/nH', 'nMg/nH', 'nSi/nH', 'nS/nH', 'nCa/nH', 
                      'nFe/nH']
    box_conversion(filename,header)    


def box_conversion(filename, header):
        
    hdf5file = filename.replace('txt','h5')

    # Read in the data
    try:
        data = np.loadtxt(filename, skiprows=2)
    except IOError:
        print 'Unable to open {0:s} to convert to HDF5'.format(filename)
        print 'Exitting'
        sys.exit()
    
    # Write data to HDF file
    try:
        df = pd.DataFrame(data, columns=header)
        df.to_hdf(hdf5file, 'data', mode='w')
    except ValueError:
        print 'Value Error with converting {0:s} in gasbox_to_hdf5'.format(filename)
        pass
    del df
    del data





def genSummaries(run,ions):

    '''
    A controlling function to generate LOS summary files.
    Uses the functions in losSummary.py

    Output file has the name:
        <galID>_a<expn>_i<incline>_<ion>_cellSummary.txt
    '''

    # Get the IDs of cells along the LOS
    probbedIDsAll = []
    for i in range(run.nlos):
        losnum = str(i+1).zfill(4)
        probbedIDsAll.append(ls.num_along_los(losnum))

    # Read in the impact parameters and azimuthal angles
    imp, phi = np.loadtxt('lines.info', skiprows=2,
                        usecols=(1,2), unpack=True)
    
    for ion in ions:

        outfile = '{0:s}_a{1:s}_i{2:d}_{3:s}_cellSummary.h5'.format(
                    run.galID,run.expn,int(run.incline),ion.name)

        print 'Compiling summary for {0:s}'.format(ion.name)
        # Write the header
        header = ['LOS','Impact','Phi','Probbed','Probbed_in_halos',
                    'In_Lines','Lines_in_halos','Significant',   
                    'Sig_in_halos']
        #f = open(outfile, 'w')
        #f.write(header)
                
        s = ('{0:d}\t{1:.3f}\t{2:.3f}\t{3:d}\t{4:d}\t{5:d}\t{6:d}\t'
                '{7:d}\t{8:d}\n')
        numcols = 9
        data = np.zeros((run.nlos,numcols))
        # Loop over the lines of sight
        for i in range(run.nlos):
            
            probbedIDs = probbedIDsAll[i]
            losnum = str(i+1).zfill(4)

            # Get the id of cells in the .lines file for this LOS 
            linesIDs = ls.num_in_lines(run.galID, ion.name, losnum)

            # Get the id of the significant cells along this LOS
            sigIDs = ls.num_significant(run.galID, run.expn, ion.name, 
                        run.incline, losnum)

            # Get the coordinates of these cells
            probbedCells, linesCells, sigCells = ls.get_coordinates(run.galID,
                run.expn, probbedIDs, linesIDs, sigIDs)

            # Determine the number of these cells that are in subhalos
            probbedinSub, linesinSub, siginSub = ls.num_in_subhalos(run.galID,
               run.expn, ion.name, run.incline, losnum, probbedCells, linesCells, 
               sigCells)

            data[i,0] = i+1
            data[i,1] = imp[i]
            data[i,2] = phi[i]
            data[i,3] = len(probbedIDs)
            data[i,4] = probbedinSub
            data[i,5] = len(linesIDs)
            data[i,6] = linesinSub
            data[i,7] = len(sigIDs)
            data[i,8] = siginSub
            # Format the output string
            #f.write(s.format(i+1, imp[i], phi[i], len(probbedIDs), 
                #probbedinSub, len(linesIDs), linesinSub, len(sigIDs),
                #siginSub))

        #f.close()
        df = pd.DataFrame(data, columns=header)
        df.to_hdr(outfile, 'data', mode='w')








