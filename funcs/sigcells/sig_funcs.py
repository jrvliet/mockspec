
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
#    f = open(linesfile)
#    redshift = float(r.readline())
#    f.close()    

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


    #################################################################
    #                                                               #
    #     Generate a noise-less spectra for the full .lines file    #
    #                                                               #
    #################################################################
#    spec.gen_spec( zabs, cellz, cellN, cellb, cellID, ion, vmax, inst
#                    transName, lamb0, fosc, gamma)
    

    # Create a los.list file containing only this LOS
#    datfile = linesfile.replace('lines', 'dat') + '\n'
#    f_los = open('los_single.list', 'w')
#    f_los.write(datfile)
#    f_los.close()
    
    # Create a Mockspec.runpars file with SNR set to zero
#    f_runpars_old = open('Mockspec.runpars')
#    f_runpars_new = open('Mockspec_0SNR.runpars', 'w')
#    l = f_runpars_old.readline().split()
#    s = l[0]+'\t\t'+l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\t'+l[7]+'\n'
#    f_runpars_new.write(s)
#    for line in f_runpars_old:
#        l = line.split()
#        l[5] = '0.'
#        s = l[0]+'\t'+l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\t'+l[7]+'\n'
#        f_runpars_new.write(s)
#    f_runpars_new.close()
#    f_runpars_old.close()

    #################################################################
    #                                                               #
    #     Generate a .lines file containing only cells with         #
    #     peculiar velocity within the velocity limits from         #
    #     the sysabs file                                           #
    #                                                               #
    #################################################################

    # New .lines file
#    newlinesfile = linesfile+'.velcut'
#    f_newlines = open(newlinesfile, 'w')
    
    # Write the first line
#    f_newlines.write('{0:.16f}\n'.format(redshift))
    
    # Loop through the original .lines file

    velz, velN, velb, velID = [], [], [], []

    velcutCount = 0
    for i in range(0,len(cell_z)):

        # Calcuate the peculiar velocity of the cell
        vpec = c*( (cell_z[i]-redshift) / (1+redshift) )
        # If the cell is inside the velocity range, write to file
        if vpec>neg_vel_limit and vpec<pos_vel_limit:
            velz.append(cell_z)
            velN.append(cell_N)
            velb.append(cell_b)
            velID.append(cell_ID)
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
    
    # Get the info from the filename
#    galID  = linesfile.split('.')[0]
#    ion    = linesfile.split('.')[1]
#    losnum = (linesfile.split('.')[2]).split('los')[1]

    # Copy the lines file for protection
#    command = 'cp '+linesfile +' '+linesfile+'.orig'
#    sp.call(command, shell=True)
#    command = 'cp '+linesfile +'.velcut '+linesfile+'.velcut.orig'
#    sp.call(command, shell=True)

    # Read in the linesfile
#    f_lines = open(linesfile+'.velcut')
#    redshift = float(f_lines.readline())
    
#    cell_z = []   # LOS redshift of each cell
#    cell_N = []   # Column density fo each cell
#    cell_b = []   # Doppler b parameter of each cell
#    cell_ID = []  # ID number of each cell
    
#    for line in f_lines:
#        l = line.split()
#        cell_z.append(float(l[0]))
#        cell_N.append(float(l[1]))
#        cell_b.append(float(l[2]))
#        cell_ID.append(int(l[3]))
#    f_lines.close()

#    if testing==1:
#        print 'In sigcells, number of velcut cells read in: ', len(cell_z)

    # Open the sysabs file
    sysabsfile  =  linesfile.replace('lines', 'sysabs')
    f_sysabs = open(sysabsfile)
    f_sysabs.readline()
    line = f_sysabs.readline()
    EW_sysabs = float(line.split()[3])
    f_sysabs.close()

    

    #################################################################
    #                                                               #
    #     Generate a noise-less spectra for the velcut .lines file  #
    #                                                               #
    #################################################################

    # Create a los.list file containing only this LOS
#    datfile = linesfile.replace('lines', 'dat') + '\n'
#    f_los = open('los_single.list', 'w')
#    f_los.write(datfile)
#    f_los.close()
    
    # Create a Mockspec.runpars file with SNR set to zero
#    f_runpars_old = open('Mockspec.runpars')
#    f_runpars_new = open('Mockspec_0SNR.runpars', 'w')
#    l = f_runpars_old.readline().split()
#    s = l[0]+'\t\t'+l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\t'+l[7]+'\n'
#    f_runpars_new.write(s)
#    for line in f_runpars_old:
#        l = line.split()
#        l[5] = '0.'
#        s = l[0]+'\t'+l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\t'+l[7]+'\n'
#        f_runpars_new.write(s)
#    f_runpars_new.close()
#    f_runpars_old.close()

    # Rename the velcut .lines to remove velcut from name, so it will be used by specsynth
#    command = 'cp '+linesfile+'.velcut '+linesfile
#    sp.call(command, shell=True)
    # Run specsynth on the velcut lines list
#    specsynth_command = codeLoc+'/funcs/mkspec/specsynth los_single.list Mockspec_0SNR.runpars'
#    sp.call(specsynth_command, shell=True)


    # Generate a noise-less spectrum for the full velcut lines
    cutLamb, cutVel, cutFlux = spec.spectrum(zabs, cutz, vutN, cutb, cutID, 
                                             ion, vmax, inst, transName, 
                                             lamb0, fosc, gamma)

    # Get the EW of this noise-less spectra
    bluewave, redwave = wavelegth(ion)
    
#    specfile = galID+'.'+ion+'.los'+losnum+'.'+bluewave+'.spec'
#    redspecfile = galID+'.'+ion+'.los'+losnum+'.'+redwave+'.spec'

    # Copy the initial quiet spectra
#    command = 'cp '+specfile+' '+specfile+'.velcutclean'
#    sp.call(command, shell=True)
#    command = 'cp '+redspecfile+' '+redspecfile+'.velcutclean'
#    sp.call(command, shell=True)
    
#    specdata = np.loadtxt(specfile)
#    wavelength = specdata[:,0]
#    velocity = specdata[:,1]
#    flux = specdata[:,2]
    
    cutEW = findEW(cutLamb, cutVel, cutFlux)
    
    # Get the goal EW
    # This is ewcut (a percent) of cutEW
    goalEW = ewcut*cutEW
    ew = cutEW

#    ewdiff = abs( (EW_sysabs - ew) / EW_sysabs )

    # Sort the cut list in order of descending N 
    rawCells = zip(cutz, cutN, cutb, cutID)
    cells = sorted(rawCells, key=itemgetter(1), reverse=True)

    # Determine which quartile to start searching

    numCells = len(cells)
    midpoint = numCells/2
    maxInd = numCells
    minInd = 0

    maxIterations = int(np.log2(numCells))
    for i in range(maxIterations):

        midpoint = int((maxInd - minInd) / 2)
        topCells = cells[:midpoint]

        # Get the ew of a spectrum using only the top 
        topz = [i[0] for i in topCells]
        topN = [i[1] for i in topCells]
        topb = [i[2] for i in topCells]
        topID = [i[3] for i in topCells]
        topLamb, topVel, topFlux = spec.gen_spec(zabs, topz, topN, topb, topID, ion,
                                                 vmax, inst, transName, lamb0, fosc,
                                                 gamma)
        topEW = findEW(topLamb, topVel, topFLux)
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
                                                   transName, lamb0, fosc, gamma)
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


    #################################################################
    #                                                               #
    #     Remove cells until the EW differs from full .lines        #
    #     value by more than ewcut                                  #
    #                                                               #
    #################################################################
    
    # Read in .lines.velcut file
#    velcutdata = np.loadtxt(linesfile+'.velcut', skiprows=1)
#    if testing==1:
#        print 'Velcutdata: ', velcutdata
#        shap = velcutdata.shape
#        print 'Velcutdata shape: ', shap
#        print 'Number of rows: ', shap[0]
#        print 'Length of shap: ', len(shap)
#        print 'Number of columns: ', shap[1]
#        print 'Index 0: ', velcutdata[1]

#    if len(cutN)==1:
        # Only one cell in file
#        velcut_z = [velcutdata[0]]
#        velcut_N = [velcutdata[1]]
#        velcut_b = [velcutdata[2]]
#        velcut_ID = [velcutdata[3]]
#    else:
#        velcut_z = list(velcutdata[:,0])
#        velcut_N = list(velcutdata[:,1])
#        velcut_b = list(velcutdata[:,2])
#        velcut_ID = list(velcutdata[:,3])
#    
#    f_log = open('sigcells.log', 'w')
#    f_log.write('Numcells \t EW \t EWdiff\n')
#    f_log.write('{0:d} \t \t{1:0.3f} \t {2:0.3f}\n'.format(len(velcut_z), ew, ewdiff))
#
#    
#    ewdiff = 0.5*ewcut
#    
#
#    
#
#
#
#    loopcount = 0
#    while ewdiff < ewcut:
#        loopcount += 1
#        
#        # Check that there is still at least one cell left
#        if len(velcut_z)>0:
#
#            # Find the cell with the lowest column denstiy
#            index = velcut_N.index(min(velcut_N))
#            
#            # Delete this index
#            del velcut_z[index]
#            del velcut_N[index]
#            del velcut_b[index]
#            del velcut_ID[index]
#
#            # Write the new .lines values to file
#            f_newlines = open(linesfile, 'w')
#            f_newlines.write('{0:.16f}\n'.format(redshift))
#            for i in range(0,len(velcut_z)):
#                s = '{0:.7f}'.format(velcut_z[i]).rjust(8)+'\t'
#                s += str(velcut_N[i]).rjust(8)+'\t'
#                s += str(velcut_b[i]).rjust(8)+'\t'
#                s += str(velcut_ID[i]).rjust(8)+'\n'
#                f_newlines.write(s)  
#            f_newlines.close()
#
#
#            # Run specsynth again
#            sp.call(specsynth_command, shell=True)
#
#            # Find the new EW 
#            specdata = np.loadtxt(specfile)
#            wavelength = specdata[:,0]
#            velocity = specdata[:,1]
#            flux = specdata[:,2]
#            
#            ew = findEW(wavelength, velocity, flux, neg_vel_limit, pos_vel_limit)
#            ewdiff = abs( (ew_velcut_lines - ew) / ew_velcut_lines)*100
#            f_log.write('{0:d} \t \t{1:0.3f} \t {2:0.3f}\n'.format(
#                        len(velcut_z), ew, ewdiff))
#        else:
#            print 'Absorption dominated by one cell'
#            singleCellCount += 1
#            ewdiff = ewcut
        
#    f_log.close()

    # Copy the lines file for protection
#    command = 'cp '+linesfile+' '+linesfile+'.final'
#    sp.call(command, shell=True)

#    return singleCellCount


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
