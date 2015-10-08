

import numpy as np


def voigt(x, y):

    """
    Computes the relative flux value in a given pixel assuming a 
    Voigt function of the Lorenzian width of the line and a 
    Gaussian broadening function; no covolution is performed

    Computes the complex probability funciton w(z)

        w(z) = exp(-z**2) * erfc(-iz) = u(x,y) + iv(x,y)
    
    where:
        z = x + iy
        u(x.y) = Voigt "function" = H(x,y)
        v(x.y) = complex component to w(z)
        x = wave number scale in Doppler width units
        y = ration of Lorentz to Doppler width

    Computation is valid for the upper half plane (y>0); 
    computational accuracy is claimed to be maximum relative error 
    < 2.0e-6 for u and <5.0e-6 for v; code is from method of 
    J. Humlicek (1979) J. Quant. Spectrosc. Radiat. Transfer _21_ 309.
    
    Accepts:
        x
        y
        
    Returns:
        u       
        v

    """
    
    m = 6

    # Load the data
    # Harmonic roots and weights of the n-point Gauss-Hermite formula
    # chosen to keep absolute error of u and v <1.0e-7
    b1 = 0.85
    t = [0.314240376,  0.947788391,  1.59768264, 2.27950708,  3.02063703,  3.8897249]
    c = [1.01172805, -0.75197147, 1.2557727e-2, 1.00220082e-2, -2.42068135e-4, 5.00848061e-7]
    s = [1.393237,  0.231152406, -0.155351466, 6.21836624e-3,  9.19082986e-5, -6.27525958e-7]

    # Intialize this pass and define the x<18.1 REGION boundary
    u = 0.0
    v = 0.0
    y1 = y + 1.5
    y2 = y1*y1
    b2 = 18.1*y + 1.65

    if y>b1 or abs(x)<b2:
        
        # Eq.(6);  REGION 1 y > (x-1.65)/18.1  for x<=18.1
        #                   y > 0.85           for x>18.1
        for i in range(0,m):

            r = x - t[i]
            d = 1.0 / (r*r+y2)
            d1 = y1*d
            d2 = r*d
    
            r = x + t[i]
            d = 1.0 / (r*r+y2)
            d3 = y1*d
            d4 = r*d

            u = u + c[i]*(d1+d3) - s[i]*(d2-d4)
            v = v + c[i]*(d2+d4) + s[i]*(d1-d3)
    else:
        
        # Eq.(11); REGION 2 y <= (x-1.65)/18.1 for x<=18.1
        #                   y <= 0.85          for x>18.1

        if abs(x)<12:
            u = np.exp(-x*x)
        y3 = y+3
        for i in range(0,m):
            r  = x - t[i]
            r2 = r*r
            d  = 1.0/(r2+y2)
            d1 = y1*d
            d2 = r*d
            u  = u + y*(c[i]*(r*d2-1.50*d1)+s[i]*y3*d2)/(r2+2.250)

            r  = x + t[i]
            r2 = r*r
            d  = 1.0/(r2+y2)
            d3 = y1*d
            d4 = r*d
            u  = u + y*(c[i]*(r*d4-1.50*d3)-s[i]*y3*d4)/(r2+2.250)
            v  = v + c[i]*(d2+d4) + s[i]*(d1-d3)

    return u, v


    

