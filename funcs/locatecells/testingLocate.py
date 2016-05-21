
import locatecells as lc
from subprocess import call

ion = 'HI'
galID = 'vela2b-28'
expn = '0.490'
sigcellsCut = 5
codeLoc = '/home/jacob/research/code/mockspec/'
inc = 90
test = 1

com='{0:s}/funcs/mklos/los7 los_single.list'.format(codeLoc)

call(com, shell=True)
print 'Starting locateSigCells\n\n'
lc.locateSigCells(galID, expn, ion, sigcellsCut, codeLoc, inc, testing=test)

 
