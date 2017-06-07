import fileinput
import numpy as np
import random
## from scipy import stats
## from numba import jit

#==============================================================================
#==============================================================================

def write_header(headstr,filenm):
    # Opens an existing file (filenm) and inserts a header at the top
    # needs fileinput to be imported
    
    for line in fileinput.input([filenm],inplace=True):
        if fileinput.isfirstline():
            print headstr
        print line,
    return 0

#==============================================================================
#==============================================================================

## @jit
def shiftarr(x,y): #Shift data to look like histogram
    # This version assumes that the x value is the center of the bin
    # instead of the left edge.
    # Assumes equal bin sizes too -- be careful!!

    # New arrays with twice the length
    newx = [None]*(len(x)*2+1)
    newy = [None]*(len(y)*2+1)

    # Actually, half of the bin width
    binwidth = np.fabs((x[0]-x[1])/2.0)

    # Start new array with first value
    newx[0] = x[0]-binwidth
    newy[0] = y[0]

    # Move over in x, but stay at same y for now
    newx[1] = x[0]+binwidth
    newy[1] = y[0]

    # counter for where you are in newx and newy
    k = 1

    # Continue with the rest of the values
    for j in range(len(y)):

        # shift x value over and y value up/down
        newx[k] = x[j]+binwidth
        newy[k] = y[j]

        k += 1

        # if you're almost to the end of your y array
        if j < (len(y)-1):
            newx[k] = x[j]+binwidth
            newy[k] = y[j+1]
            k += 1

    # Stay at same y value, shift over in x
    newx[-1] = x[j]+binwidth
    newy[-1] = 0

    # return your new arrays
    return newx,newy

#==============================================================================
#==============================================================================

## @jit
def bindata(binsample,binmin,binmax,binsize):
    # CWC's bindata-linear with normalization by number of data points.
    # no errors here - those are calculated by bootstrap/jackknife

    mindat = binmin
    maxdat = binmax+binsize
    nbins = int((maxdat - mindat) / binsize)

    # initialize new arrays
    bincounts = []
    binslo,binshi = [],[]
    bincenter = []

    # Loop through each bin
    for i in range(0,nbins):

        # initial bins
        dumbins = 0.0
        binlo = mindat + i*binsize
        binhi = binlo + binsize

        # determine bin counts given bins
        for j in range(len(binsample)):
            if binsample[j] >= binlo and binsample[j] < binhi:
                dumbins += 1.0
            elif binsample[j] == nbins and binsample[j] == binhi:
                dumbins += 1.0

        # normalize by total number of points
        dumbins = dumbins/len(binsample)

        # Fill new arrays with values calculated in this loop
        bincounts.append(dumbins)
        binslo.append(binlo)
        binshi.append(binhi)
        bincenter.append(binlo + (0.5*binsize))

    return bincenter,bincounts

#==============================================================================
#==============================================================================

## @jit
def deltaVelocity(vdat):
    # takes velocity data (vdat) and calculates the difference between
    # every possible pair
    
    dVi = []
    for i in range(len(vdat)):
        for j in range(i+1,len(vdat)):
            dVi.append(abs(vdat[i]-vdat[j]))

    return dVi

#==============================================================================
#==============================================================================


def bootstrap(galdict,zgals,Vcirc):
    # galdict is the dictionary of the galaxies/regions you want to
    # bootstrap, e.g., blue galaxies. Determining the sample of
    # galaxies must be done outside of this function

    # galdict = {galkey1:([vdats1]),galkey2:([vdats2])...}
    
    # zgals and Vcirc are just used to determine the bin sizes

    # Full sample velocity separations + array of keys for random draw
    galkeys = [] 
    dVi = []

    # loop through galaxy keys in the dictionary
    for galaxy in galdict:
        galkeys.append(galaxy)

        # clear vdat array from last galaxy
        vdat = []

        # for each key, grab the array of velocities, put them all in
        # new array. So vdat contains pixel velocities for A SINGLE
        # galaxy (absorber) in your sample/dictionary.
        for velocity in galdict[galaxy]:
            vdat.append(velocity)

        # Calculate velocity separations for a single galaxy
        # (absorber) in sample. Append velocity splittings to dVi for
        # binning
        ## dVi = deltaVelocity(vdat) # old line
        dumdeltav = deltaVelocity(vdat)
        for val in dumdeltav:
            dVi.append(val)

    # Max dVi for full sample
    maxdVi = max(dVi)

    # calculate bin sizes. If not worrying about zgal or Vcirc, just
    # use binsize=10 (or some other number) and you don't need to pass
    # zgals or Vcirc into this function - change def
    # bootstrap(galdict,zgals,Vcirc) above to def bootstrap(galdict)
    if zgals[0] == 0 and Vcirc[0] == 0: binsize = 10.0
    elif zgals[0] != 0 and Vcirc[0] == 0: binsize = 33.33

    elif zgals[0] == 0 and Vcirc[0] != 0: binsize = 0.1
    elif zgals[0] != 0 and Vcirc[0] != 0: binsize = 0.3333


    # Bin up the full data sample - This is your TPCF!
    bincenter,bincounts = bindata(dVi,0.0,maxdVi,binsize)


    # Array to keep track of the binned bootstrap realizations
    #
    # Will be a 2D array where each line is a different boostrap
    # realization and each column is a different TPCF velocity
    # separation bin
    bincountruns = []

    # Loop through galaxies/absorbers to do bootstrap analysis
    for i in range(10): # I do 1000
        print i

        # array for pixel velocities for this bootstrap realization
        dumdVi = []

        # randomly draw absorbers with replacement
        for j in range(len(galdict)):
            #randgal is the index for galkeys
            randgal = random.randrange(0,len(galdict),1)

            # clear newvdat array from last galaxy
            newvdat = []

            # for each key, grab the array of velocities, put them all
            # in new array. only one absorber at a time
            for velocity in galdict[galkeys[randgal]]:
                newvdat.append(velocity)

            # Calculate velocity separations for this absorber. total
            # dumdVi will represent one bootstrap realization.
            ## dumdVi = deltaVelocity(newvdat) # old line
            dumdeltav = deltaVelocity(newvdat)
            for val in dumdeltav:
                dumdVi.append(val)

        # bin the data and keep track of values in each bin
        dumcenter,dumcounts = bindata(dumdVi,0.0,maxdVi,binsize)

        # Append bootstrap realization binned data to array
        bincountruns.append(dumcounts)

    # calculate uncertainties in each bincountruns bin.
    # axis=0 does mean and std along columns
    meanbins = np.mean(bincountruns,axis=0)
    stdbins = np.std(bincountruns,axis=0)
    
    return bincenter,bincounts,dVi,meanbins,stdbins


#==============================================================================
#==============================================================================

## @jit
def updownerrors(binmean,binstd):
    # Convert standard devations to up and down errors for TPCF
    # plotting

    binstddn,binstdup = [],[]

    for i in range(len(binmean)):
        binstddn.append(binmean[i] - binstd[i])
        binstdup.append(binmean[i] + binstd[i])

    return binstddn,binstdup


#==============================================================================
#==============================================================================

## @jit
def makedictionary(dirloc,zgal,Vcirc):
    # Make dictionary for the subsample using directory locations
    # (dirloc), zgal, and Vcirc

    cvel = 2.99792458*10.**5. #km/s

    kindir = '/home/matrix2/nnielsen/Kinematics/'
    #kindir = '../'

    # if using an equivalent width sensitivity cut, read in the cut value
    ewcut = np.genfromtxt(kindir+'Tables/TPCF/region_ewcut.dat')

    # just the directory where all of the dirloc data live
    kiniso = kindir+'Kin/iso/'

    # initialize the subsample dictionary
    ssdict = {}

    # loop through the directory locations
    for i in range(len(dirloc)):
        vdat = []

        # get the absorber redshift (zabs)
        zabs = np.genfromtxt(kiniso+dirloc[i]+'/zabs.dat',usecols=0)

        # get the equivalent widths of each region to compare to the
        # EW sensitivity cut
        reg,regEWr = np.genfromtxt(kiniso+dirloc[i]+'/MgII2796.ews',
                                   skip_header=1,usecols=(0,4),unpack=True)

        # get the region numbers and velocity ranges
        regnum,regvm,regvp = np.genfromtxt(kiniso+dirloc[i]+'/MgII2796.ewreg',
                                   skip_header=1,usecols=(0,5,6),unpack=True)

        # Get the spectrum - pixel velocities
        pixelvel = np.genfromtxt(kiniso+dirloc[i]+'/MgII2796',usecols=1)

        # Check each region EW to see if it is greater than the
        # sensitivity cut
        validreg = []
        for j in range(len(np.atleast_1d(regEWr))):

            # stupid zero dimensional arrays require this nonsense
            if regEWr.shape != ():
                dumregEWr = regEWr[j]
                dumreg    = reg[j]
            elif regEWr.shape == ():
                dumregEWr = float(regEWr[()])
                dumreg    = float(reg[()])

            # If greater, then keep!
            if dumregEWr >= ewcut:
                validreg.append(dumreg)

        # for each valid region, get the velocity ranges
        for j in range(len(validreg)):
            
            # stupid zero dimensional arrays require this nonsense
            if len(np.atleast_1d(regvm)) == 1:
                velminus = regvm
                velplus  = regvp
            elif len(np.atleast_1d(regvp)) > 1:
                velminus = regvm[validreg[j]-1]
                velplus  = regvp[validreg[j]-1]

            # loop through each pixel and check to see if it is in the
            # valid region
            for pix in pixelvel:
                if pix > velminus and pix < velplus:

                    # Keep this line for your stuff, but drop the if
                    # statement
                    if Vcirc[i] == 0 and zgal[i] == 0:
                        vdat.append(pix)

                    # Normalize pixel velocity by the circular
                    # velocity if Vcirc given
                    elif Vcirc[i] != 0 and zgal[i] == 0:
                        vdat.append(pix/Vcirc[i])

                    # If zgal given, shift velocities to zgal instead
                    # of zabs. You don't need this statement for your
                    # stuff either
                    elif zgal[i] != 0:

                        zshift = (zabs - zgal[i])/(1. + zgal[i])
                        zshift = zshift*cvel

                        vdum = abs(pix + zshift)

                        # Don't normalize by Vcirc
                        if Vcirc[i] == 0:
                            vdat.append(vdum)

                        # Normalize by Vcirc
                        elif Vcirc[i] != 0:
                            vdat.append(vdum/Vcirc[i])

        # Since I'm randomly drawing regions for bootstrap instead
        # of galaxies (or absorbers), give the dictionary keys a
        # name with dirloc and the region number
        ## dummyreg,dum2 =  str(validreg[j]).split('.') # old line
        ## galreg = dirloc[i]+'_'+dummyreg # old line

        # Make dictionary now of -- {dirloc+region_number:[pixel vels]}
        ## ssdict[galreg] = vdat # old line

        # Make dictionary of absorbers (galaxies), not regions!
        ssdict[dirloc[i]] = vdat

    return ssdict
