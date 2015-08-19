
import numpy as np
import sys
import subprocess as sp
from ew import findEW 

def velcut(linesfile):

    """
    Cuts the cells out of the .lines file that do not fall withing
    the velocity window as found in the sysabs file
    """
#    import numpy as np
#    import sys
#    from  math import sqrt, pow, log10
#    import os.path
#    import os
#    import subprocess as sp
    
    # Define constants
    c = 3.0e5   # Speed of light in km/s

    # Get the info from the filename
    galID  = linesfile.split('.')[0]
    ion    = linesfile.split('.')[1]
    losnum = (linesfile.split('.')[2]).split('los')
    
    # Read in the linesfile
    f_lines = open(linesfile)
    redshift = float(f_lines.readline())

    cell_z = []   # LOS redshift of each cell
    cell_N = []   # Column density fo each cell
    cell_b = []   # Doppler b parameter of each cell
    cell_ID = []  # ID number of each cell

    for line in f_lines:
        l = line.split()
        cell_z.append(float(l[0]))
        cell_N.append(float(l[1]))
        cell_b.append(float(l[2]))
        cell_ID.append(int(float(l[3])))
    f_lines.close()
    
    # Open the sysabs file
    sysabsfile  =  linesfile.replace('lines', 'sysabs')
    f_sysabs = open(sysabsfile)
    f_sysabs.readline()
    line = f_sysabs.readline()
    neg_vel_limit = float(line.split()[1])
    pos_vel_limit = float(line.split()[2])
    EW_sysabs = float(line.split()[3])
    f_sysabs.close()

    #################################################################
    #                                                               #
    #     Generate a noise-less spectra for the full .lines file    #
    #                                                               #
    #################################################################

    # Create a los.list file containing only this LOS
    datfile = linesfile.replace('lines', 'dat') + '\n'
    f_los = open('los_single.list', 'w')
    f_los.write(datfile)
    f_los.close()
    
    # Create a Mockspec.runpars file with SNR set to zero
    f_runpars_old = open('Mockspec.runpars')
    f_runpars_new = open('Mockspec_0SNR.runpars', 'w')
    l = f_runpars_old.readline().split()
    s = l[0]+'\t\t'+l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\t'+l[7]+'\n'
    f_runpars_new.write(s)
    for line in f_runpars_old:
        l = line.split()
        l[5] = '0.'
        s = l[0]+'\t'+l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\t'+l[7]+'\n'
        f_runpars_new.write(s)
    f_runpars_new.close()
    f_runpars_old.close()

    #################################################################
    #                                                               #
    #     Generate a ,lines file containing only cells with         #
    #     peculiar velocity within the velocity limits from         #
    #     the sysabs file                                           #
    #                                                               #
    #################################################################

    # New .lines file
    newlinesfile = linesfile+'.velcut'
    f_newlines = open(newlinesfile, 'w')
    
    # Write the first line
    f_newlines.write('{0:.16f}\n'.format(redshift))
    
    # Loop through the original .lines file
    for i in range(0,len(cell_z)):

        # Calcuate the peculiar velocity of the cell
        vpec = c*( (cell_z[i]-redshift) / (1+redshift) )
        
        # If the cell is inside the velocity range, write to file
        if vpec>neg_vel_limit and vpec<pos_vel_limit:
            s = '{0:.7f}'.format(cell_z[i]).rjust(8)+'\t'
            s += str(cell_N[i]).rjust(8)+'\t'
            s += str(cell_b[i]).rjust(8)+'\t'
            s += str(cell_ID[i]).rjust(8)+'\n'
            f_newlines.write(s)

    f_newlines.close()


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def sigcells(linesfile, ewcut):
    # Description:
    #  Determines which cells are the significant contributers 
    #  to the EW measurement

    
    # Get the info from the filename
    galID  = linesfile.split('.')[0]
    ion    = linesfile.split('.')[1]
    losnum = (linesfile.split('.')[2]).split('los')[1]

    # Copy the lines file for protection
#    command = 'cp '+linesfile +' '+linesfile+'.orig'
#    sp.call(command, shell=True)
#    command = 'cp '+linesfile +'.velcut '+linesfile+'.velcut.orig'
#    sp.call(command, shell=True)

    # Read in the linesfile
    f_lines = open(linesfile+'.velcut')
    redshift = float(f_lines.readline())
    
    cell_z = []   # LOS redshift of each cell
    cell_N = []   # Column density fo each cell
    cell_b = []   # Doppler b parameter of each cell
    cell_ID = []  # ID number of each cell
    
    for line in f_lines:
        l = line.split()
        cell_z.append(float(l[0]))
        cell_N.append(float(l[1]))
        cell_b.append(float(l[2]))
        cell_ID.append(int(l[3]))
    f_lines.close()

    # Open the sysabs file
    sysabsfile  =  linesfile.replace('lines', 'sysabs')
    f_sysabs = open(sysabsfile)
    f_sysabs.readline()
    line = f_sysabs.readline()
    neg_vel_limit = float(line.split()[1])
    pos_vel_limit = float(line.split()[2])
    EW_sysabs = float(line.split()[3])
    f_sysabs.close()



    #################################################################
    #                                                               #
    #     Generate a noise-less spectra for the velcut .lines file  #
    #                                                               #
    #################################################################

    # Create a los.list file containing only this LOS
    datfile = linesfile.replace('lines', 'dat') + '\n'
    f_los = open('los_single.list', 'w')
    f_los.write(datfile)
    f_los.close()
    
    # Create a Mockspec.runpars file with SNR set to zero
    f_runpars_old = open('Mockspec.runpars')
    f_runpars_new = open('Mockspec_0SNR.runpars', 'w')
    l = f_runpars_old.readline().split()
    s = l[0]+'\t\t'+l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\t'+l[7]+'\n'
    f_runpars_new.write(s)
    for line in f_runpars_old:
        l = line.split()
        l[5] = '0.'
        s = l[0]+'\t'+l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\t'+l[7]+'\n'
        f_runpars_new.write(s)
    f_runpars_new.close()
    f_runpars_old.close()

    # Rename the velcut .lines to remove velcut from name, so it will be used by specsynth
    command = 'cp '+linesfile+'.velcut '+linesfile
    sp.call(command, shell=True)
    # Run specsynth on the velcut lines list
    specsynth_command = '/home/matrix2/cwc/Projects/Mockspec/Codes/mkspec/specsynth los_single.list Mockspec_0SNR.runpars'
    sp.call(specsynth_command, shell=True)

    # Get the EW of this noise-less spectra
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

    specfile = galID+'.'+ion+'.los'+losnum+'.'+bluewave+'.spec'
    redspecfile = galID+'.'+ion+'.los'+losnum+'.'+redwave+'.spec'

    # Copy the initial quiet spectra
    command = 'cp '+specfile+' '+specfile+'.velcutclean'
    sp.call(command, shell=True)
    command = 'cp '+redspecfile+' '+redspecfile+'.velcutclean'
    sp.call(command, shell=True)
    
    specdata = np.loadtxt(specfile)
    wavelength = specdata[:,0]
    velocity = specdata[:,1]
    flux = specdata[:,2]
    
    ew = findEW(wavelength, velocity, flux, neg_vel_limit, pos_vel_limit)
    ewdiff = abs( (EW_sysabs - ew) / EW_sysabs )
    ew_velcut_lines = ew

    #################################################################
    #                                                               #
    #     Remove cells until the EW differs from full .lines        #
    #     value by more than ewcut                                  #
    #                                                               #
    #################################################################

    # Read in .lines.velcut file
    velcutdata = np.loadtxt(linesfile+'.velcut', skiprows=1)
    velcut_z = list(velcutdata[:,0])
    velcut_N = list(velcutdata[:,1])
    velcut_b = list(velcutdata[:,2])
    velcut_ID = list(velcutdata[:,3])
    
    f_log = open('sigcells.log', 'w')
    f_log.write('Numcells \t EW \t EWdiff\n')
    f_log.write('{0:d} \t \t{1:0.3f} \t {2:0.3f}\n'.format(len(velcut_z), ew, ewdiff))

    
    ewdiff = 0.5*ewcut
    

    while ewdiff < ewcut:
        
        # Find the cell with the lowest column denstiy
        index = velcut_N.index(min(velcut_N))
        
        # Delete this index
        del velcut_z[index]
        del velcut_N[index]
        del velcut_b[index]
        del velcut_ID[index]

        # Write the new .lines values to file
        f_newlines = open(linesfile, 'w')
        f_newlines.write('{0:.16f}\n'.format(redshift))
        for i in range(0,len(velcut_z)):
            s = '{0:.7f}'.format(velcut_z[i]).rjust(8)+'\t'
            s += str(velcut_N[i]).rjust(8)+'\t'
            s += str(velcut_b[i]).rjust(8)+'\t'
            s += str(velcut_ID[i]).rjust(8)+'\n'
            f_newlines.write(s)  
        f_newlines.close()


        # Run specsynth again
        sp.call(specsynth_command, shell=True)

        # Find the new EW 
        specdata = np.loadtxt(specfile)
        wavelength = specdata[:,0]
        velocity = specdata[:,1]
        flux = specdata[:,2]
        
        ew = findEW(wavelength, velocity, flux, neg_vel_limit, pos_vel_limit)
        ewdiff = abs( (ew_velcut_lines - ew) / ew_velcut_lines)
        f_log.write('{0:d} \t \t{1:0.3f} \t {2:0.3f}\n'.format(
                    len(velcut_z), ew, ewdiff))
        
    f_log.close()

    # Copy the lines file for protection
    command = 'cp '+linesfile+' '+linesfile+'.final'
    sp.call(command, shell=True)




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
