

# Code to test the spectrum code
import matplotlib.pyplot as plt
import spectrum as sp
import files as fi
import sys

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


diff = []
filebase = './files/vela29v2_1.CIV.los'
for i in range(1,1000):
    losnum = '{0:04}'.format(i)
    filename = filebase+losnum+'.lines'
    
    # Read in test los
    f = open(filename)
    zabs = float(f.readline())

    zcell, Ncell, bcell, cellID = [], [], [], []
    for line in f:
        l = line.split()
        zcell.append(float(l[0]))
        Ncell.append(float(l[1]))
        bcell.append(float(l[2]))
        cellID.append(float(l[3]))
    f.close()

    lamb, vel, flux = sp.spectrum(zabs, zcell, Ncell, bcell, cellID, 
                                  ion, vmax, inst, transName, 
                                  lamb0, fosc, gamma)
    

    # Read in control
    specFile = filename.replace('lines', 'CIV1548.spec')
    f = open(specFile)
    clamb, cvel, cflux = [], [], []
    for line in f:
        l = line.split()
        clamb.append(float(l[0]))
        cvel.append(float(l[1]))
        cflux.append(float(l[2]))
    f.close()

    # Compare results
    for i in range(0,len(flux)):
        diff.append(abs(flux[i] - cflux[i]))
    
# Print results
print 'Mean difference: {0:.6f}'.format(np.mean(diff))
print 'Min difference:  {0:.6f}'.format(np.min(diff))
print 'Max difference:  {0:.6f}'.format(np.max(diff))
print 'Std difference:  {0:.6f}'.format(np.std(diff))


###############3
# Plot


f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

specfile = 'test.spec'
f = open(specfile)
vel, wave, flux = [], [], []
for line in f:
    vel.append(float(line.split()[1]))
    wave.append(float(line.split()[0]))
    flux.append(float(line.split()[2]))
ax1.step(wave, flux, 'k', label='raw')
ax1.plot(wave, flux, 'kx')
f.close()

testflux = flux
testwave = wave


specfile = 'test.convolve'
f = open(specfile)
vel, wave, flux = [], [], []
for line in f:
    vel.append(float(line.split()[1]))
    wave.append(float(line.split()[0]))
    flux.append(float(line.split()[2]))
ax2.step(wave, flux, 'g', label='convolve')
ax2.plot(wave, flux, 'gx')
f.close()



specfile = './files/vela29v2_1.CIV.los0001.CIV1548.spec'
f = open(specfile)
vel, wave, flux = [], [], []
for line in f:
    vel.append(float(line.split()[1]))
    wave.append(float(line.split()[0]))
    flux.append(float(line.split()[2]))
ax3.step(wave, flux, 'r', label='control')
ax3.plot(wave, flux, 'rx')
f.close()

#plt.xlim([3093, 3101])
ax1.legend(frameon=False, loc='lower right')
ax2.legend(frameon=False, loc='lower right')
ax3.legend(frameon=False, loc='lower right')
ax1.set_ylabel('Flux')
ax2.set_ylabel('Flux')
ax3.set_ylabel('Flux')
ax3.set_xlabel('Wavelength')
plt.savefig('test.pdf')
plt.savefig('test.jpg')

print 'Length of test wave: {0:d}'.format(len(testwave))
print 'Length of control wave: {0:d}'.format(len(wave))

print 'Test wave step = {0:f}'.format(testwave[1]-testwave[0])
print 'Control step = {0:f}'.format(wave[1]-wave[0])



diff = 0.0
total = 0.0
for i in range(0,len(flux)):

    diff += testflux[i] - flux[i]
    total += flux[i]

err = diff / total 
print 'Percent error: {0:%}'.format(err)

plt.cla()
plt.clf()


f = open('./files/fresponse.dat')
fres = []
for line in f:
    fres.append(float(line))
f.close()

f = open('response.dat')
res = []
for line in f:
    res.append(float(line))
f.close()

plt.plot(fres, 'rx', label='control')
plt.plot(res, 'gx', label='test')
plt.legend(frameon=False, loc='lower right')
plt.savefig('isf.pdf')

print 'Length of response: ', len(fres), len(res)
print 'Sum of response: ', sum(fres), sum(res)
