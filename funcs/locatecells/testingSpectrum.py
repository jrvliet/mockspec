

# Code to test the spectrum code
import matplotlib.pyplot as plt
import spectrum as sp
import files as fi

ion = 'CIV'
inst = 'COSNUV'
vmax = 1000

# Get transition info
k,j,transName,lamb0,fosc,gamma,mamu,abund,ip = fi.get_transitions(ion)

print '\nFrom get_transitions:'
print '\tk         = {0:d}'.format(k)
print '\tj         = {0:d}'.format(j)
print '\ttransName = {0:s}'.format(transName)
print '\tlamb0     = {0:f}'.format(lamb0)
print '\tfosc      = {0:f}'.format(fosc)
print '\tgamma     = {0:f}'.format(gamma)
print '\tmamu      = {0:f}'.format(mamu)
print '\tabund     = {0:f}'.format(abund)
print '\tip        = {0:f}'.format(ip)

# Read in test los
filename = 'vela29v2_1.CIV.los0001.lines'
f = open(filename)
zabs = float(f.readline())

zcell, Ncell, bcell, cellID = [], [], [], []
for line in f:
    l = line.split()
    zcell.append(float(l[0]))
    Ncell.append(float(l[1]))
    bcell.append(float(l[2]))
    cellID.append(float(l[3]))

sp.spectrum(zabs, zcell, Ncell, bcell, cellID, ion, vmax, inst, 
            transName, lamb0, fosc, gamma)


specfile = 'test.spec'
f = open(specfile)
vel, wave, flux = [], [], []
for line in f:
    vel.append(float(line.split()[1]))
    wave.append(float(line.split()[0]))
    flux.append(float(line.split()[2]))

plt.step(wave, flux, 'k', label='test')
plt.plot(wave, flux, 'kx')
f.close()
print 'Test Wave Step: {0:f} starting at {1:f}'.format(wave[1]-wave[0], wave[0]) 

testflux = flux


specfile = 'vela29v2_1.CIV.los0001.CIV1548.spec'
f = open(specfile)
vel, wave, flux = [], [], []
for line in f:
    vel.append(float(line.split()[1]))
    wave.append(float(line.split()[0]))
    flux.append(float(line.split()[2]))

plt.step(wave, flux, 'r', label='control')
plt.plot(wave, flux, 'rx')
f.close()
print 'Control Wave Step: {0:f} starting at {1:f}'.format(wave[1]-wave[0], wave[0]) 

plt.xlim([3095, 3097.5])
plt.legend(frameon=False, loc='lower left')
plt.savefig('test.pdf')


diff = 0.0
total = 0.0
for i in range(0,len(flux)):

    diff += testflux[i] - flux[i]
    total += flux[i]

err = diff / total 
print 'Percent error: {0:%}'.format(err)
