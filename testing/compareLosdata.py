import sys

errorcount = 0
trueloc = './good/'
losfiles = 'losdata.list'
files = open(losfiles)
for line in files:

    newf = open(line.strip())
    oldf = open(trueloc+line.strip())
    
    for i in range(2):
        newf.readline()
        oldf.readline()

    for newline, oldline in zip(newf, oldf):

        newl = newline.split()
        oldl = oldline.split()
        for i in range(0,len(newl)):

            if float(newl[i])!=float(oldl[i]):
    
                errorcount += 1
print errorcount    
