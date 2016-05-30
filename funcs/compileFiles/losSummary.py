
'''
Creates a file that summarizes the results of each LOS in terms
of the cells

Contains:
    - number of cells probbd by the LOS
    - number of cells included in the .lines file for each ion
    - number of cells that contribute significantly to the absorption
    - number of cells located within subhalos
    - number of cluster
'''

from math import sqrt
from numpy import loadtxt


def num_along_los(losnum):

    '''
    Determines the number of cells probed by the LOS.
    Determiened by the length of the cellID.dat file output by cellfinder
    There is no header
    '''


    filename = 'los{0:s}.cellID.dat'.format(losnum)
    fileloc = './cellIDs/{1:s}'.format(ion,filename)
    
    # Get the number of cells along the LOS
    ids = []
    with open(fileloc, 'r') as f:
        for line in f:
            ids.append(int(line.strip()))

    return ids



def num_in_lines(galID, ion, losnum):

    '''
    Determines the number of cells included in the .lines file.
    To be included in the .lines file, the cell must have a column
    density contribution of at least 10^9
    '''

    headerSize = 1

    filename = './{0:s}/{1:s}.{0:s}.los{2:s}.lines'.format(ion,galID,losnum)

    # Get the number of cells it the file
    ids = []
    with open(filename, 'r') as f:
        f.readline()
        for line in f:
            ids.append(int(line.strip().split()[3]))

    return ids



def num_significant(galID, expn, ion, inc, losnum):
    
    '''
    Determines the number of cells along the LOS that contribute significantly
    to the absorption. Determined by running through the abs_cells.dat file
    and counting those lines that have that losnum
    '''
    
    filename = '{0:s}/{1:s}.{2:s}.{0:s}.i{3:s}.abs_cells.dat'.format(ion,galID,expn,ion,inc)

    linecount = 0
    found = 0
    
    ids = []
    with open(filename, 'r') as f:
        # Read past the header
        f.readline()
        for line in f:
            los = int(line.split()[0])
            if los == int(losnum):
                ids.append(int(line.split()[2]))
                if found == 0:
                    found = 1
            else:
                if found == 1:
                    break

    return ids

def get_coordinates(galID, expn, probbedIDs, linesIDs, sigIDs):

    '''
    Gets the x, y, z coordinates of each cell based on the id
    '''
   
    gasfile = '../{0:s}_GZa{1:s}.txt'.format(galID, expn)
    xgas, ygas, zgas = loadtxt(gasfile, skiprows=2, usecols=(1,2,3), unpack=True)
     
    # Loop through probbed cells
    xprob, yprob, zprob = [], [], []
    for cell in probbedIDs:
        index = cell-1
        xprob.append(xgas[index])
        yprob.append(ygas[index])
        zprob.append(zgas[index])
    
    # Loop through lines cells
    xlines, ylines, zlines = [], [], []
    for cell in probbedIDs:
        index = cell-1
        xlines.append(xgas[index])
        ylines.append(ygas[index])
        zlines.append(zgas[index])
    
    # Loop through sig cells
    xsig, ysig, zsig = [], [], []
    for cell in probbedIDs:
        index = cell-1
        xsig.append(xgas[index])
        ysig.append(ygas[index])
        zsig.append(zgas[index])
    
    probbedCells = (xprob, yprob, zprob)
    linesCells = (xlines, ylines, zlines)
    sigCells = (xsig, ysig, zsig)

    return probbedCells, linesCells, sigCells

def num_in_subhalos(galID, expn, ion, inc, losnum, 
                    probbedCells, linesCells, sigCells):

    '''
    Determines the number of cells along the LOS that reside within
    the halo of a satellite galaxy. Only include the 10 largest
    satellite galaxies
    '''

    numHalos = 10

    # Read in the subhalo information
    halofile = '../input_{0:s}.txt'.format(expn)

    xhalo, yhalo, zhalo, rhalo = [], [], [], []
    with open(halofile, 'r') as f:
        # Read past the main halo
        l = f.readline().split()
        xhost = float(l[0])
        yhost = float(l[1])
        zhost = float(l[2])

        for i in range(numHalos):      
            l = f.readline().split()
            xhalo.append( (float(l[0])-xhost)*1000.0  )
            yhalo.append( (float(l[1])-yhost)*1000.0  )
            zhalo.append( (float(l[2])-zhost)*1000.0  )
            rhalo.append( (float(l[4]) )

    
    # Run through the probbed cells
    insubCount = 0
    for cell in probbedCells:
        x, y, z = cell

        insub = 0
        # Loop through the subhalos
        for i in range(numHalos):
            
            # Get the distance between the cell and this halo
            r = sqrt( (x-xhalo[i])**2 + (y-yhalo[i])**2 + (z-zhalo[i])**2 ) 
            if r<rhalo:
                insub = 1
    
        insubCount += insub

    probbedinSub = insubCount

    # Run through the lines cells
    insubCount = 0
    for cell in linesCells:
        x, y, z = cell

        insub = 0
        # Loop through the subhalos
        for i in range(numHalos):
            
            # Get the distance between the cell and this halo
            r = sqrt( (x-xhalo[i])**2 + (y-yhalo[i])**2 + (z-zhalo[i])**2 ) 
            if r<rhalo:
                insub = 1
    
        insubCount += insub

    linesinSub = insubCount

    # Run through the significant cells
    insubCount = 0
    for cell in sigCells:
        x, y, z = cell

        insub = 0
        # Loop through the subhalos
        for i in range(numHalos):
            
            # Get the distance between the cell and this halo
            r = sqrt( (x-xhalo[i])**2 + (y-yhalo[i])**2 + (z-zhalo[i])**2 ) 
            if r<rhalo:
                insub = 1
    
        insubCount += insub

    siginSub = insubCount

    
    return probbedinSub, linesinSub, siginSub













