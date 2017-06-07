
class runProps(object):

    '''
    A class to describe the parameters of the run of mockspec. 
    Basically holds everything read in from mockspec.config
    in read_control_file except ion information
    '''
    
    def __init__ (self):
        
        # General properties
        self.galID = 'galID'
        self.expn = '1.000'
        self.nlos = 1000
        self.maximpact = 1.5
        self.incline = 0
        self.ewcut = 0.0
        self.snr = 30
        self.ncores = 1
        self.rootLoc = ''
        self.sigcellsCut = 5.
        
        # Flags
        self.runRates = 1
        self.runGenLOS = 1
        self.runCellfinder = 1
        self.runIdcells = 1
        self.runLos7 = 1
        self.runSpecsynth = 1
        self.runSysanal = 1
        self.runCullabs = 1
        self.runLocateCells = 1
        self.runSummaries = 1
        self.runPlotting = 1
        
        
class ionProps(object):

    '''
    Class to describe an ion being studied. Atributtes include the 
    ion name, the metallicity enhancement, and the instrument
    '''
    
    def __init__ (self):
        
        self.name = 'HI'
        self.xh = 0.
        self.instrument = 'COSNUV'
    
