# Functions for manipulating files for generating
# spectra


def config_spec(instrument):

    """
    Reads Mockspec.instruments and get the properties for
    the desired instrument
    
    Accepts:
        instrument (name of instrument)

    Returns: 
        R_isf (resolution)
        presel (pixels per resolution element)
        RN (read noise)
    """

    filename = 'Mockspec.instruments'
    with open(filename) as f:
        
        # Read past header:
        f.readline()

        # Loop through file and locate the instrument name
        for line in f:
            


    
    
    
