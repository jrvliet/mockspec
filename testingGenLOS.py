

from genLOS import *

galID = 'vela29v2'
gasfile = 'vela29v2_GZa0.500.txt'
summaryLoc = '/lustre/projects/p089_swin/jvander/mockspec/summaries/'
expn = '0.500'
inc = 75
nLOS = 2
maximpact = 1.5
ncores = 1
doPlot =1 

for i in range(0,90,1):
    genLines(galID, gasfile, summaryLoc, expn, i, nLOS, maximpact, ncores, doPlot)




