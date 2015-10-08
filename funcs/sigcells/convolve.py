
import numpy as np
from scipy import signal as sg
from scipy import interpolate as si



def convolve(wrkflux, response, nfft, ncondat, ndata, resfac, mode='same'):

    # Interpolate wrkflux onto convdata. Only take up ncondat
    # of the array
    x = np.linspace(0, 10, num=len(wrkflux))
    y = wrkflux
    f = si.interp1d(x, y, kind='cubic')

    xnew = np.linspace(0, 10, num=ncondat)
    ynew = f(xnew)
    
    highResFlux = np.ones(nfft)
    for i in range(0,len(ynew)):
        highResFlux[i] = ynew[i]


    # Pad the beginning
    x = np.ones(len(highResFlux + 100))
    for i in range(100,len(x)):
        x[i] = highResFlux[i]

    flux = sg.fftconvolve(highResFlux, response, mode=mode)

    # Pick off every resfac element and stuff into idx element of flux
    for i in range(100, ndata):
        k = (i-1)*resfac
        wrkflux[i] = abs(flux[int(k)])


    return wrkflux

