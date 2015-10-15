

#def findEW(wavelength, velocity, flux, neglimit, poslimit):
def findEW(wavelength, velocity, flux):
    # Get step in wavelength
    step = wavelength[1] - wavelength[0]

    ew = 0.0
    
    # Loop through wavelengths
    for i in range(0,len(wavelength)):
        
#        if velocity[i]>neglimit and velocity[i]<poslimit:
            
        midflux = (flux[i]+flux[i+1]) / 2.0
        ew += (1.0-midflux)*step


    return ew
