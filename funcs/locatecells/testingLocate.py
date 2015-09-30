
import locatecells as lc
from subprocess import call

ion = 'HI'
galID = 'vela20'
expn = '0.440'
sigcellsCut = 5
codeLoc = '/lustre/projects/p089_swin/jvander/mockspec/'
test = 1

com='/lustre/projects/p089_swin/jvander/mockspec/funcs/mklos/los7 los_single.list'

call(com, shell=True)
print 'Starting locateSigCells\n\n'
lc.locateSigCells(galID, expn, ion, sigcellsCut, codeLoc, testing=test)

 
