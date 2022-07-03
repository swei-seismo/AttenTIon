#!/opt/antelope/5.4/bin/python
## PROGRAM TO GET P/S ARRIVALS FOR ANTELOPE DATABASE
## Written by S. Wei, MARCH 2015

import os,sys
sys.path.append(os.environ['ANTELOPE'] +'/data/python')
from antelope.datascope import *
import globaldb as g

workdir='/P/weisq/attentomo/testdb'
if not os.path.isdir(workdir):
    os.mkdir(workdir)

##s=raw_input('Usage: mkgsfl.py subdb# orid')
##subdb=s.split()[0]
##orid=int(s.split()[1])

for subdb in ['1']:
    evtlst='/P/weisq/attentomo/input/eventid_sub'+subdb+'.lst'
    oridlst=[int(line.split()[0]) for line in open(evtlst).readlines()[1:]]
    for orid in oridlst:
        print(subdb,orid)

        outpfln=workdir+'/arr%s_%d.txt' % (subdb,orid)
        outsfln=workdir+'/arr%s_%d.txt' % (subdb,orid)
        outpfile=open(outpfln,'w')
        outsfile=open(outsfln,'w')

        g.dbread(subdb)

        dbev=g.dbass.subset('orid==%d' % orid)
        dbev=dbev.sort('sta')
        num =dbev.query('dbRECORD_COUNT')
##        print(num)

        for ii in range(num):
            dbev.record=ii
            (sta,time,phs)=dbev.getv('sta','time','iphase')
##            print(sta,phs)
##            depth=max(depth,0)
##            print(sta,delta,depth,phs)
            if not sta=='' or phs=='':
                if phs=='P':
                    outpfile.write('%s %s %f\n' % (sta,phs,time))
                elif phs=='S':
                    outsfile.write('%s %s %f\n' % (sta,phs,time))
                else:
                    print('Wrong phase: %s' % phs)
                    continue

        outpfile.close()
        outsfile.close()
