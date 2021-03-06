# Functions for creating the model absorption for spectrum.py

import numpy as np
#from scipy import constants as sc
from astropy import constants as const
import voigt as vp

def set_atomic(lamb0, f0, gamma0):

    c = const.c.cgs.value
    e = const.e.esu.value
    m_e = const.m_e.cgs.value    

    con1Const = ( np.sqrt(np.pi)*e*e ) / ( m_e*c*c )
    con2Const = 4 * np.pi * c 
    con1 = 1e5 * f0 * lamb0*lamb0 * con1Const
    con2 = 1e-8 * gamma0 * lamb0*lamb0 / con2Const

    return con1, con2




def init_spectrum(waveMin, waveMax, dwave):
    """
    Initilizes the spectrum to proper pixelization 
    and unity value
    
    Accepts:
        waveMin (minimum wavelength in the spectrum to make)
        waveMax (maximum wavelength in the spectrum to make)
        dwave (pixel sampling rate)

    Returns:
        ndata (number of data points)
        lamb (array of wavelegth values)
        wrkflux (array of flux values, all set to 1)

    """

    # Given the wavelength range and the pixel sampling rate, 
    # determine the number of pixels "ndata" in the spectrum
    ndata = int( (waveMax - waveMin)/dwave ) + 1

    # Fill the wavelength pixel values and initialize the 
    # spectrum as unitu continuum in array wrkflux

    lamb = np.zeros(ndata)
    wrkflux = np.ones(ndata)
    for i in range(0,ndata):
        lamb[i] = waveMin + (i)*dwave

    return ndata, lamb, wrkflux




def do_abs_lines(lamb0, zabs, nline, zline, bline, con1, con2, lamb, wrkflux):


    """
    Uses Voigt profile formalism and backgorund source
    radiative transfer (no source function) to compute
    the absorption lines
    
    Some helpful definitions:
    b1  = column density
    b2  = rest frame central wavelength
    b3  = rest frame Doppler width
    y   = natural broadening in terms of doppler units
    x   = wavelength difference in doppler units
    w0  = rest frame wavelegnth for transition
    w   = rest frame wavelength fro observed wavelength
    tau = optical depth


    Accepts:
        lamb0 (central wavelength of transition)
        zabs  (redshift of absorption)
        nline (array of cell column densities)
        zline (array of cell redshifts)
        bline (array of cell Doppler b values)
        con1  (constatn 2 form set_atomic)
        con2  (constant 1 from set_atomic)
        lamb (array of wavelengths)
        wrkflux (array of flux values)

    Returns:
        wrkflux (array of flux values)

    """

    ckms = const.c.to('km/s').value

    for i in range(0,len(nline)):

        # Initalize Voight profile parameters
        # These must be rest from quantities even though 
        # the line may be redshifted
        w0 = lamb0
        b1 = nline[i]
        b2 = w0 * (1.0+zline[i]) / (1.0+zabs)
        b3 = w0 * abs(bline[i])/ckms
        y  = con2 / b3
        tcon = con1 * (b1/b3)
#        print w0, b1, b2, b3, y, tcon
        # Loop over the pixels and perform the radiative transfer
        for j in range(0,len(wrkflux)):
            w = lamb[j] / (1.0+zabs)
            x = (w-b2)/b3
            u, v = vp.voigt(x,y)
            wrkflux[j] = wrkflux[j] * np.exp(-tcon*u)

    return wrkflux
    


#def do_lyman_limit():
#def drop_lines():

# The following functions were originally found in instrument.f
def instrument(ndata, wcen, flag, R_fac, dwave):
    """
    Sets up the intrumental spread function for convolving the
    original model spectrum

    Given the instrumental profile sigma in velocity units;
    the ISF profiles are loaded in wrap-around order for the 
    convolution; zero spatial information is in the first index
    
    flag = 0 ; call for setting up convolution
    flag = 1 ; =0 + communicates intrumental parameters

    Accepts:
        ndata (number of data points)
        wcen (center something)
        flag (a flag as defined above)

    Returns:
        response (instrumental ISF)

    """
    
    # Following parameters are set in files.f in the original verison
    # but there are no places where they ever get modified,
    # so in the interest of avoiding long funciton calls, 
    # set them here
    slit = 1.0
    resfac = 3.0
    conwindo = 3


    np2 = 17
    ckms = const.c.to('km/s').value
    pwrsof2 = [16,32,64,128,256,512,1024,2048,4096,8192,
               16384,32768,65536,131072,262144,524288,1048576]

    # Compute the instrumental resolution in velocity units [km/s]
    # sigma in km/s given by FWHM/2.35
    # profile = sigma of the Gaussian ISF in velocity [km/s]
    if (R_fac==0.0):
        profile = slit
    else:
        profile = slit*ckms / (2.35*R_fac)

    # dv is the velocity sampling of the pixels km/s/pixel...
    # the number of pixels per resolution element = profile/dv
    dv = dwave/wcen * ckms
    pixpres = 2.35 * profile/(dv*slit)
    hdv = dv/resfac
    nresponse = 2*int(conwindo*profile/hdv) + 1
    
    # Stuff the response function
    response = np.zeros(nresponse)
    response[0] = phi(0, profile)
    norm = response[0]

    for i in range(0,nresponse):
        xdv = i*hdv
        response[i] = phi(xdv, profile)
        norm += response[i]

    # Currently the response is one sided. Fold it over
    isf = []
    for i in range(0,len(response)):
        isf.append(response[-1*i])
    for i in range(0,len(response)):
        isf.append(response[i])
    isf[0] = isf[-1]
    norm = norm*2 
    # For the convolution integral, the response function needs
    # to be normalized or the flux is not conserved
    # Unit width is used becuase data are discretized by pixel indices 
    # in the convolution
    norm = sum(response)
    for i in range(0,nresponse):
        response[i] /= norm
    isf = [i/sum(isf) for i in isf]    


    # Compute the length of the convolution functions
    ncondat = int(resfac)*(ndata-1) - 1
    nfft = ncondat + (nresponse-1)/2 + 1
    for i in range(0,np2):
        if nfft<pwrsof2[i]:
            nfft = pwrsof2[i]
            break

    f = open('response.dat', 'w')
    for i in isf:
        f.write('{0:.10f}\n'.format(i))
    f.close()
        

    return isf, nfft, ncondat


def phi(dvw, width):

    """
    Given the sigma in velocity units, this routine computes
    the relative Gaussian value of the instrumental profile
    Called by routine instrument iteratively
    """

    # dvw is the value at which the instrumental response is to be evaluated
    z = dvw/width
    val = np.exp(-1.0*(z*z)/2.0)

    return val








def calc_velocity(lamb, zabs, lamb0):

    """
    Create the velocity array for outputing to the .spec file
    """

    ckms = const.c.to('km/s').value
    velocity = []
    for i in range(0,len(lamb)):
        wave = lamb[i]
        w0z = (1.0+zabs) * lamb0
        vel = ckms*(wave-w0z)/w0z
        velocity.append(vel)

    return velocity



def convolve(flux, ):

    """
    The spectrum is convolved with the instrumental
    spread function to make the resulting model
    spectrum. 

    Accepts: 


    Returns:
        wrkflux: Flux after convolution
    """







