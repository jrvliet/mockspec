#!/usr/bin/python


fin = open('D9o2.wrong')
fout = open('D9o2.dat', 'w')

i = 0
for line in fin:
    i += 1
    if i==1 or i==2:
        fout.write(line)
    else:
        l = line.split()
        for j in range(4,len(l)):
            val = float(l[j])/10.0
            l[j] = '{0:0.6f}'.format(val)

#        s = '{0:.3f} \t {1:0.5f} \t {2:0.4e} \t {3:.5f} \t {4:.8f} \t {5:.8f} \t {6:.8f} \t {7:.8f} \t {8:.8f} \t {9:.8f} \t {10:.8f} \t {11:.8f} \t {12:.8f} \t {13:.5f} \t {14:.5f} \t {15:.5f} \n'.format(l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12], l[13], l[14], l[15])
        s = ''
        for j in range(0,len(l)):
            s += l[j]+'\t'
        s += '\n'
        fout.write(s)
fin.close()
fout.close
