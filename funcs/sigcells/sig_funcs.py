
import numpy as np
import sys
import subprocess as sp
from ew import findEW 
import spectrum as spec
from operator import itemgetter

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

    # Generate a noise-less spectrum for the full velcut lines
    cutLamb, cutVel, cutFlux = spec.gen_spec(zabs, cutz, cutN, cutb, cutID, 
                                             ion, vmax, inst, transName, 
                                             lamb0, fosc, gamma)

    # Get the EW of this noise-less spectra
    bluewave, redwave = wavelength(ion)
    
    cutEW = findEW(cutLamb, cutVel, cutFlux)
    
    # Get the goal EW
    # This is ewcut (a percent) of cutEW
    goalEW = ewcut*cutEW
    ew = cutEW

    # Sort the cut list in order of descending N 
    rawCells = zip(cutz, cutN, cutb, cutID)
    cells = sorted(rawCells, key=itemgetter(1), reverse=True)

    # Determine which quartile to start searching

    numCells = len(cells)
    midpoint = numCells/2
    maxInd = numCells
    minInd = 0
    print numCells

    maxIterations = int(np.log2(numCells))
    for i in range(maxIterations):

        midpoint = int((maxInd - minInd) / 2)
        topCells = cells[:midpoint]

        print len(topCells)
        # Get the ew of a spectrum using only the top 
        topz = [i[0] for i in topCells]
        topN = [i[1] for i in topCells]
        topb = [i[2] for i in topCells]
        topID = [i[3] for i in topCells]
        topLamb, topVel, topFlux = spec.gen_spec(zabs, topz, topN, topb, topID, ion,
                                                 vmax, inst, transName, lamb0, fosc,
                                                 gamma)
        topEW = findEW(topLamb, topVel, topFlux)
        if topEW>goalEW:
            minInd = midpoint
        else:
            maxInd = midpoint

        
    
    # Verify this cut point
    roughCells = cell[:midpoint]
    roughz = [i[0] for i in roughCells]
    roughN = [i[1] for i in roughCells]
    roughb = [i[2] for i in roughCells]
    roughID = [i[3] for i in roughCells]
    roughLamb, roughVel, roughFlux = spec.gen_spec(zabs, roughz, roughN, roughb,
                                                   roughID, ion, vmas, inst, 
                                                   transName, lamb1, fosc, gamma)
    roughEW1 = findEW(roughlamb, roughVel, roughFlux)

    roughCells = cell[:midpoint+1]
    roughz = [i[0] for i in roughCells]
    roughN = [i[1] for i in roughCells]
    roughb = [i[2] for i in roughCells]
    roughID = [i[3] for i in roughCells]
    roughLamb, roughVel, roughFlux = spec.gen_spec(zabs, roughz, roughN, roughb,
                                                   roughID, ion, vmas, inst, 
                                                   transName, lamb0, fosc, gamma)
    roughEW2 = findEW(roughlamb, roughVel, roughFlux)

    # roughEW1 should be below the goalEW
    # roughEW2 shoudl be above the goalEW
    if roughEW1<goalEW and roughEW2>goalEW:
        return roughz, roughN, roughb, roughID
    else:
        raise ValueError        




def nT(filename, norm):

    
    data = np.loadtxt(filename, skiprows=1)
    
    # Do the normalized histogram
    if norm==1:
        galID = filename.split('.')[0]
        e1 = filename.split('.')[1]
        e2 = filename.split('.')[2]
        expn = e1+'.'+e2
    
        fullbox_name = '../'+galID+'_GZa'+expn+'.txt'
        
        fullbox = np.loadtxt(fullbox_name, skiprows=2)
        
        nfull = np.log10(fullbox[:,7])
        Tfull = np.log10(fullbox[:,8])

        print max(nfull), min(nfull)
        print max(Tfull), min(Tfull)
    


    n = data[:,7]
    T = data[:,8]

    
    # Bin the data
    numbins = 50
    H, xedges, yedges = np.histogram2d( n, T, bins=numbins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    
    # Rotate and filp the histogram
    H = np.rot90(H)
    H = np.flipud(H)

    # Mask the bins where the count is zero
    Hmasked = np.ma.masked_where(H==0,H)
    
    # Take the log of the count
    Hmasked = np.log10(Hmasked)

    # Plot and save the figure
    fig = plt.figure()
    plt.pcolormesh(xedges,yedges,Hmasked)
    plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
    plt.ylabel(' $\log$ (T) [K] ')
    plt.xlim([-8, 1])
    plt.ylim([2,8])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('$\log$ (Counts)')

    outname = filename.replace('dat','phase.pdf')
    plt.savefig(outname)


    print H
    if norm==1:

        print 'Norm'
        Hf, xedgesf, yedgesf = np.histogram2d( nfull, Tfull, bins=numbins)
        extentf = [xedgesf[0], xedgesf[-1], yedgesf[0], yedgesf[-1]]

        Hf = np.rot90(Hf)
        Hf = np.flipud(Hf)

        print np.max(Hf)
        print Hf
        print type(Hf)
        print Hf.shape
        Hnorm = H
        for i in range(0,len(Hf[:,0])):
            for j in range(0,len(Hf[0,:])):
                if Hf[i,j]>0:
                    Hnorm[i,j] = H[i,j] / Hf[i,j]
                else:
                    Hnorm[i,j] = 0
        
        Hnorm = np.ma.masked_where(Hnorm==0,Hnorm)

        print Hnorm
        plt.cla()
        plt.clf()
        fig = plt.figure()
        plt.pcolormesh(xedges,yedges,Hnorm)
        plt.xlabel(' $\log (n_{H})$ [cm$^{-3}$] ')
        plt.ylabel(' $\log$ (T) [K] ')
        plt.xlim([-8, 1])
        plt.ylim([2,8])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Counts')
        plt.clim(0,1)
        outname = filename.replace('dat','normphase.pdf')
        plt.savefig(outname)





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



 
