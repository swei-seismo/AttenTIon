## PROGRAM TO ADD NOISE TO SYNTHETIC DATA OF t*

import os
from numpy import random


s=raw_input('Please input infile and outfile: ')
infl=s.split()[0]
outfl=s.split()[1]
phase=s.split()[2]
phase=phase.upper()

## READ INPUT FILE
if phase=='P':
    stn =[line.split()[0] for line in open(infl).readlines()]
    lat =[float(line.split()[1]) for line in open(infl).readlines()]
    lon =[float(line.split()[2]) for line in open(infl).readlines()]
    dep =[float(line.split()[3]) for line in open(infl).readlines()]
    tst =[float(line.split()[4]) for line in open(infl).readlines()]
    ter =[float(line.split()[5]) for line in open(infl).readlines()]
    tsd =[float(line.split()[6]) for line in open(infl).readlines()]
    att =[float(line.split()[7]) for line in open(infl).readlines()]
elif phase=='S' or phase=='SP' or phase=='PS':
    stn =[line.split()[0] for line in open(infl).readlines()]
    lat =[float(line.split()[1]) for line in open(infl).readlines()]
    lon =[float(line.split()[2]) for line in open(infl).readlines()]
    dep =[float(line.split()[3]) for line in open(infl).readlines()]
    stst=[float(line.split()[4]) for line in open(infl).readlines()]
    ster=[float(line.split()[5]) for line in open(infl).readlines()]
    ptst=[float(line.split()[6]) for line in open(infl).readlines()]
    pter=[float(line.split()[7]) for line in open(infl).readlines()]
    att =[float(line.split()[8]) for line in open(infl).readlines()]
    QpQs=[float(line.split()[9]) for line in open(infl).readlines()]
else:
    exit('Wrong phase: '+phase)

## ADD NOISE AND OUTPUT
outlst=open(outfl,'w')
for i in range(len(stn)):
    if phase=='P':
        newtst=tst[i]+random.normal(0,ter[i])
        if newtst<=0:
            newtst=tst[i]-tst[i]*random.rand()
        outlst.write('%5s%12.4f%12.4f%12.4f%12.6f%5.2f%5.2f%12.6f\n'
                 %(stn[i],lat[i],lon[i],dep[i],newtst,tsd[i],ter[i],att[i]))
    elif phase=='S':
        newtst=stst[i]+random.normal(0,ster[i])
        if newtst<=0:
            newtst=stst[i]-stst[i]*random.rand()
        outlst.write('%5s%12.4f%12.4f%12.4f%12.6f%5.2f%12.6f%5.2f%12.6f%7.2f\n'
            %(stn[i],lat[i],lon[i],dep[i],newtst,ster[i],ptst[i],pter[i],
              att[i],QpQs[i]))
    elif phase=='SP' or phase=='PS':
        newstst=stst[i]+random.normal(0,ster[i])
        if newstst<=0:
            newstst=stst[i]-stst[i]*random.rand()
        newptst=ptst[i]+random.normal(0,pter[i])
        if newptst<=0:
            newptst=ptst[i]-ptst[i]*random.rand()
        outlst.write('%5s%12.4f%12.4f%12.4f%12.6f%5.2f%12.6f%5.2f%12.6f%7.2f\n'
            %(stn[i],lat[i],lon[i],dep[i],newstst,ster[i],newptst,pter[i],
              att[i],QpQs[i]))
outlst.close()
