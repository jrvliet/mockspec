
import numpy as np
import sys
import subprocess as sp
from ew import findEW 
import spectrum as spec
from operator import itemgetter
import matplotlib.pyplot as plt

def velcut(cellz, cellN, cellb, cellID, linesfile, redshift, testing=0):

    """
    Cuts the cells out of the .lines file that do not fall withing
    the velocity window as found in the sysabs file
    """
    # Define constants
    c = 3.0e5   # Speed of light in km/s

    # Get the redshift of the absorption
    # This is the first line of the the lines file

    # Open the sysabs file
    sysabsfile  =  linesfile.replace('lines', 'sysabs')
    f = open(sysabsfile)
    f.readline()
    line = f.readline()
    neg_vel_limit = float(line.split()[1])
    pos_vel_limit = float(line.split()[2])
    f.close()

    if testing==1:
        print '\t\tFrom sysabs:'
        print '\t\t\tNeg Vel Limt: {0:f}'.format(neg_vel_limit)
        print '\t\t\tPos_vel_limi: {0:f}'.format(pos_vel_limit)
        print '\t\t\tEW:           {0:f}'.format(EW_sysabs)


    velz, velN, velb, velID = [], [], [], []

    velcutCount = 0
    for i in range(0,len(cellz)):

        # Calcuate the peculiar velocity of the cell
        vpec = c*( (cellz[i]-redshift) / (1+redshift) )
        # If the cell is inside the velocity range, write to file
        if vpec>neg_vel_limit and vpec<pos_vel_limit:
            velz.append(cellz[i])
            velN.append(cellN[i])
            velb.append(cellb[i])
            velID.append(cellID[i])
            velcutCount += 1
     
    if testing==1:
        print '\t\tAfter velcut, number of cells: ', velcutCount

    return velz, velN, velb, velID

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def significant_cells(zabs, cutz, cutN, cutb, cutID, ewcut, codeLoc, 
                      transProps, testing=0):
    # Description:
    #  Determines which cells are the significant contributers 
    #  to the EW measurement

    # Unpack transProps
    ion, vmax, inst, transName, lamb0, fosc, gamma = transProps

    singleCellCount = 0       # Counts number of LOS dominated by a single cell
#    f = open('cutN_1.log', 'w')
#    for n in cutN:
#        f.write('{0:f}\n'.format(n))
#    f.close()
    # Generate a noise-less spectrum for the full velcut lines
    cutLamb, cutVel, cutFlux = spec.gen_spec(zabs, cutz, cutN, cutb, cutID, 
                                             ion, vmax, inst, transName, 
                                             lamb0, fosc, gamma)
    plt.plot(cutVel, cutFlux, label='Vel Cut')

#    f = open('cutN_2.log', 'w')
#    for n in cutN:
#        f.write('{0:f}\n'.format(n))
#    f.close()
    # Get the EW of this noise-less spectra
    bluewave, redwave = wavelength(ion)
    
    cutEW = findEW(cutLamb, cutVel, cutFlux)
    
    # Get the goal EW
    # This is ewcut (a percent) of cutEW
    ewcut = ewcut / 100.0
    goalEW = ewcut*cutEW
    ew = cutEW
#    print 'ewcut = {0:f}'.format(ewcut)
#    print 'cutEW = {0:f}'.format(cutEW)
#    print 'goalEW = {0:f}'.format(goalEW)

    # Sort the cut list in order of descending N 
    rawCells = zip(cutz, cutN, cutb, cutID)
    cells = sorted(rawCells, key=itemgetter(1), reverse=True)

    
#    f = open('cutN_3.log', 'w')
#    for c, n, b, id in cells:
#        f.write('{0:f}\n'.format(n))
#    f.close()
    # Determine which quartile to start searching

    numCells = len(cells)
    midpoint = numCells/2
    maxInd = numCells
    minInd = 0
    print numCells

    maxIterations = int(np.log2(numCells))

#    flog = open('ew.log', 'w')
#    flog.write('EW = {0:f}\tNum Cells = {1:d}\n'.format(ew, numCells))
#    for c, n, b, id in cells:
#        flog.write('{0:f}\n'.format(n))
#    flog.write('\n')
 
#    print 'Max Iterations = ', maxIterations
    for i in range(maxIterations):

        midpoint = int((maxInd - minInd) / 2)
        topCells = cells[:midpoint]

#        print '\nMin Index = {0:d} \t Max Index = {1:d} \t Midpoint = {2:d} \t Length of topCells: {3:d}'.format(minInd, maxInd, midpoint, len(topCells))
        # Get the ew of a spectrum using only the top 
        topz = [j[0] for j in topCells]
        topN = [j[1] for j in topCells]
        topb = [j[2] for j in topCells]
        topID = [j[3] for j in topCells]
        topLamb, topVel, topFlux = spec.gen_spec(zabs, topz, topN, topb, topID, ion,
                                                 vmax, inst, transName, lamb0, fosc,
                                                 gamma)
        topEW = findEW(topLamb, topVel, topFlux)

#        print 'i = ',i
#        plt.plot(topVel, topFlux, label='Cut {0:d}'.format(i+1))

#        print 'EW = {0:f}\tGoal = {1:f}\tNum Cells = {2:d}'.format(topEW, goalEW, len(topN))
        

#        flog.write('EW = {0:f}\tNum Cells = {1:d}\n'.format(topEW, len(topN)))
#        for c, n, b, id in topCells:
#            flog.write('{0:f}\t{1:f}\n'.format(n, c))
#        flog.write('\n')

        if topEW>goalEW:
            maxInd = midpoint
        else:
            minInd = midpoint

#    flog.close()
#    plt.legend(frameon=False)
#    plt.savefig('ew.pdf')
    # Verify this cut point
    roughCells = cells[:midpoint]
    roughz = [i[0] for i in roughCells]
    roughN = [i[1] for i in roughCells]
    roughb = [i[2] for i in roughCells]
    roughID = [i[3] for i in roughCells]
    roughLamb, roughVel, roughFlux = spec.gen_spec(zabs, roughz, roughN, roughb,
                                                   roughID, ion, vmax, inst, 
                                                   transName, lamb0, fosc, gamma)
    roughEW1 = findEW(roughLamb, roughVel, roughFlux)
#    print 'Length of rough 1 = {0:d}'.format(len(roughN))
    numRoughCells = len(roughN)
    roughCells = cells[:midpoint+1]
    roughz = [i[0] for i in roughCells]
    roughN = [i[1] for i in roughCells]
    roughb = [i[2] for i in roughCells]
    roughID = [i[3] for i in roughCells]
    roughLamb, roughVel, roughFlux = spec.gen_spec(zabs, roughz, roughN, roughb,
                                                   roughID, ion, vmax, inst, 
                                                   transName, lamb0, fosc, gamma)
    roughEW2 = findEW(roughLamb, roughVel, roughFlux)
#    print 'Length of rough 2 = {0:d}'.format(len(roughN))

    # roughEW1 should be below the goalEW
    # roughEW2 shoudl be above the goalEW
#    print 'roughEW1 = {0:f}\nroughEW2 = {1:f}\ngoalEW = {2:f}'.format(roughEW1,
#                                                                      roughEW2,
#                                                                      goalEW)
#    print numRoughCells
    if (roughEW1<goalEW and roughEW2>goalEW) or numRoughCells==1:
        return roughz, roughN, roughb, roughID
    else:
        raise ValueError        





def wavelength(ion):
    if ion == 'HI':
        bluewave = 'Lya'
        redwave = ion+'1026'        
    elif ion == 'MgII':
        bluewave = ion+'2796'
        redwave = ion+'2803'
    elif ion == 'CIV':
        bluewave = ion+'1548'
        redwave = ion+'1551'
    elif ion == 'OVI':
        bluewave = ion+'1032'
        redwave = ion+'1038'
    else:
        print 'Unkown ion'
        sys.exit()

    return bluewave, redwave




def get_transition_properties( ion, codeLoc):

    f = open( codeLoc+'/data/Mockspec.transitions')
    f.readline()
    found = 0
    for line in f:
        l = line.split()
        flag = int(l[0])
        species = l[4]
        if ion==species and flag==1 and found==0:
            found = 1
            element = l[1]
            j = l[3]
            transName = l[4]
            lamb0 = float(l[6])
            fosc = float(l[7])
            gamma = float(l[8])
    f.close()

    f = open('Mockspec.runpars')
    f.readline()
    for line in f:
        l = line.split()
        if element in l[0] and j==l[1]:
            inst = l[3]
            vmax = float(l[6])
    f.close()

    return (ion, vmax, inst, transName, lamb0, fosc, gamma)         



 
