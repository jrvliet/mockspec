
""" Identifies which cells are the significant source of absorpiton
    """


import numpy as np
import subprocess as sp
import ew 
import os.path as op
import spectrum as spec
import sig_funcs as sf


def sigCells(galID, expn, ion, ewcut, codeLoc, testing=0):

    redshift = 1.0/float(expn) - 1.0
    singleCount = 0     # Counts number of LOS dominated by a single cell 

    # Get the properties of the transition
    transProps = sf.get_transition_properties(ion, codeLoc)
    
    
    # Read in the galaxy's box
    boxfile = galID+'_GZa'+expn+'.'+ion+'.txt'
    box = np.loadtxt(boxfile, skiprows=2)
    if testing==1:
        print 'Box read in'

    # Read in the LOS info from lines.info
    losInfo = np.loadtxt('lines.info',skiprows=2)

    # Open the output file
    out_file = galID+'.'+expn+'.'+ion+'.abs_cells.dat'
    f_out = open(out_file, 'w')

    # Write a header to the output file
    header = 'LOS \t Imp Param      Cell ID   Redshift        logN'\
            '    Doppler b    Galactocentric d      log nH     log T'\
             '     Cell Size     SNII mass frac      SNIa mass frac'\
             '       alpha_Zmet\n'
    f_out.write(header)



    # Get a list of LOS that have a sysabs file associated with it
    sysabsLosnum = []
    for i in range(0,1000):
        losnum = str(i).zfill(4)
        filename = galID+'.'+ion+'.los'+losnum+'.sysabs'
        # Check to see if the file exists
        if op.isfile(filename):
           sysabsLosnum.append(losnum)

    if testing==1:
         print 'Sysabs files aggregated'

    for i in range(0,len(sysabsLosnum)):

        losnum = sysabsLosnum[i]
        num = int(losnum)

        if testing==1:
            print 'LOS num: ', losnum

        linesfile = galID+'.'+ion+'.los'+losnum+'.lines'
    
        # Read in the lines file
        f = open(linesfile)
        zabs = float(f.readline())

        cellz, cellN, cellb, cellID = [], [], [], []
        for line in f:
            l = line.split()
            cellz.append(float(l[0]))
            cellN.append(float(l[1]))
            cellb.append(float(l[2]))
            cellID.append(float(l[3]))
        f.close()
        numcells = len(cellz)

        # Check that there are cells along the LOS
        if numcells>0:

            # Get the impact parameter of this LOS 
            imp = losInfo[num-1, 1]

            # Perform the velocity cut
            if testing==1:
                print '\t Performing velocity cut'
            cutz, cutN, cutb, cutID = sf.velcut(cellz, cellN, cellb, cellID, 
                                                linesfile, redshift, 
                                                testing=testing)

            # Find the significant cells
            if testing==1:
                print '\t Finding significant cells'
            sigz, sigN, sigb, sigID = sf.significant_cells(zabs, cutz, cutN, cutb, 
                                                           cutID, ewcut, codeLoc,
                                                           transProps, testing=testing)

            # Get the properties of the cells
            for j in range(0,len(sigz)):

                # Get the cell's properties from the boxfile
                index = sigID[j] - 1
                cellSize = box[index, 0]
                x = box[index, 0]
                y = box[index, 0]
                z = box[index, 0]
                density = np.log10(box[index, 7])
                temperature = np.log10(box[index, 8])
                snII = box[index, 9]
                snIa = box[index, 10]
                alphaZ = box[index, 16]

                # Calculate the galactocentric distance
                r = np.sqrt( x*x + y*y + z*z)

                # Write all to the output file
                s = ('{0:d}'.format(num).ljust(7) + 
                    '{0:.3f}'.format(imp).rjust(7) +
                    '{0:d}'.format(cellID).rjust(16) + 
                    '{0:-.7f}'.format(redshift).rjust(14) + 
                    '{0:.3f}'.format(column).rjust(10) + 
                    '{0:.3f}'.format(doppler).rjust(13) + 
                    '{0:.5e}'.format(r).rjust(20) + 
                    '{0:.4f}'.format(density).rjust(12) + 
                    '{0:.4f}'.format(temperature).rjust(10) + 
                    '    {0:.4f}'.format(cell_size).rjust(14) + 
                    '{0:.4e}'.format(snII).rjust(19) + 
                    '{0:.4e}'.format(snIa).rjust(20) + 
                    '{0:.4    e}'.format(alphaZ).rjust(17) + '\n' )
                f_out.write(s)


    f_out.close()

