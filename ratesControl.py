

def setupRatesControl(gasfile, expn, ion_list, requiredLoc)

    if not os.path.exists('rates.inp'):
        command = 'cp '+requiredLoc+'rates.inp .'
        sp.call(command, shell=True)

    # Alter the rates.inp file
    numions = len(ion_list)
    defaultratesinp = open('rates.inp')
    ratesinp = open('rates.inp.tmp', 'w')
    for i in range(0,3):
        line = defaultratesinp.readline()
        ratesinp.write(line)
    # Replace the gas file name
    line = defaultratesinp.readline().replace('MW9_GZ932.a1.001.txt', gasfile)
    ratesinp.write(line)
    # Replace the expansion parameter
    line = defaultratesinp.readline().replace('1.001', expn)
    ratesinp.write(line)
    # Replace the number of output cubes
    line = defaultratesinp.readline().replace('4', str(numions))
    ratesinp.write(line)
    # Copy the rest of the file
    for line in defaultratesinp:
        ratesinp.write(line)

    defaultratesinp.close()
    ratesinp.close()

    command = 'rm rates.inp'
    sp.call(command, shell=True)
    command = 'mv rates.inp.tmp rates.inp'
    sp.call(command, shell=True)


def setupRatesOutputs():

    if not os.path.exists('rates.outfiles'):
        command = 'cp '+requiredLoc+'rates.outfiles .'
        sp.call(command, shell=True)
    
    # Alter the rates.outfiles file
    defoutfile = open('rates.outfiles')
    outfiles = open('rates.outfiles.tmp', 'w')
    for ion in ion_list:
        ionbox = galID+'_GZa'+expn+'.'+ion+'.txt'
        element, Z, excitation = get_transition_info(ion)
        line = '{0:<27s} {1:>2s} {2:>2s}\n'.format(ionbox, Z, excitation)
        outfiles.write(line)
    outfiles.close()
    defoutfile.close()
    command = 'rm rates.outfiles'
    sp.call(command, shell=True)
    command = 'mv rates.outfiles.tmp rates.outfiles'
    sp.call(command, shell=True)
   

def rates():
    
    print 'Generating ion boxes...'
     

