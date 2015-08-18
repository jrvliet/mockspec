import sys
import numpy as np
errorcount = 0
trueloc = './good/'
losfiles = 'losdata.list'
files = open(losfiles)
count = 0
lines = 0
smallcount = 0
per = []
ind = []
bad = open('badcells.out', 'w')
f = open('badfiles.out', 'w')
for line in files:

    newf = open(line.strip())
    oldf = open(trueloc+line.strip())
    
    wrote = 0
    for i in range(2):
        header = newf.readline().split()
        oldf.readline()
        

    for newline, oldline in zip(newf, oldf):
        cellwrote = 0
        lines += 1
        newl = newline.split()
        oldl = oldline.split()
        for i in range(0,len(newl)):
            count += 1
            n = float(newl[i])
            o = float(oldl[i])
            if n!=o:
                errorcount += 1
                if abs(n-o)<=0.011:
                    smallcount += 1
                percent = (n-o)/o
                per.append(percent)
                ind.append(i)
                cellnum = newl[-1]
                if cellwrote==0:
                    bad.write(cellnum+'\n')
                    cellwrote = 1
                if wrote==0:
                    f.write(line)
                    wrote=1
                print 'Filename: {0:s}  Line Num: {1:d}   Line index: {2:d}   Percent: {3:.3%}   Old: {4:.3f}    New: {5:.3f}  Diff: {6:f}'.format(line.strip(), lines, i, percent, o, n, n-o)
#    print len(ind)
#    break
f.close()
ind = set(ind)

bad.close()

print 'Errors:    ', errorcount    
print 'Small err: ', smallcount
print 'Cells:  ', count
print 'Lines:  ', lines
print 'Max Percent: {0:.3%}'.format(max(per))
print 'Max Percent: {0:.3%}'.format(np.max(per))
print 'Min Percent: {0:.3%}'.format(np.min(per))
print 'Mean Percent: {0:.3%}'.format(np.mean(per))

print 'Unique Indicies: ', ind
for i in range(len(header )):
    print i, header[i]

