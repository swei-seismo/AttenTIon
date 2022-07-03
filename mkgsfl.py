#!/opt/antelope/5.4/bin/python
## PROGRAM TO CREATE GEOMETRIC SPREADING FILES FOR P AND S
## Written by S. Wei, MARCH 2015

import os,sys
sys.path.append(os.environ['ANTELOPE'] +'/data/python')
from antelope.datascope import *
import numpy as np
import globaldb as g

workdir='/P/weisq/attentomo/GS'
momcmd='/P/weisq/attentomo/program/src/momcalc'
if not os.path.isdir(workdir):
    os.mkdir(workdir)

##s=raw_input('Usage: mkgsfl.py subdb# orid')
##subdb=s.split()[0]
##orid=int(s.split()[1])

##for subdb in ['1','2','3','4','11','22','33','44']:
for subdb in ['3']:
##for subdb in ['22']:
    evtlst='/P/weisq/attentomo/input/eventid_sub'+subdb+'.lst'
    oridlst=[int(line.split()[0]) for line in open(evtlst).readlines()[1:]]
    for orid in oridlst:
##    for orid in [51879]:
##    for orid in [55961]:
        print(subdb,orid)

        outpfln=workdir+'/pgsfile%s_%d.txt' % (subdb,orid)
        outsfln=workdir+'/sgsfile%s_%d.txt' % (subdb,orid)
        outpfile=open(outpfln,'w')
        outsfile=open(outsfln,'w')

        g.dbread(subdb)

        dbev=g.dbass.subset('orid==%d' % orid)
        dbev=dbev.sort('sta')
        num =dbev.query('dbRECORD_COUNT')
##        print(num)

        inimomp=np.array([1.0e-12,5.0e-13,1.0e-13])
        inimoms=np.array([1.0e-12,5.0e-13,1.0e-13])
        for ii in range(num):
            dbev.record=ii
            (sta,delta,depth,phs)=dbev.getv('sta','delta','depth','iphase')
            depth=max(depth,0)
##            print(sta,delta,depth,phs)
            if phs=='P':
                cmd='%s -9. %f %f P' %(momcmd,delta,depth)
            elif phs=='S':
                cmd='%s -9. %f %f SH' %(momcmd,delta,depth)
            else:
                print('Wrong phase: %s' % phs)
                continue
            osout=os.popen(cmd).read()
##            print(cmd)
##            print(osout)
            moment=float(osout.split()[0])
##            print(moment)
            if phs=='P':
                if moment!=0:
                    inimomp[0]=inimomp[1]
                    inimomp[1]=inimomp[2]
                    inimomp[2]=moment
                if moment==0:
                    moment=np.mean(inimomp[np.nonzero(inimomp)])
            else:
                if moment!=0:
                    inimoms[0]=inimoms[1]
                    inimoms[1]=inimoms[2]
                    inimoms[2]=moment
                if moment==0:
                    moment=np.mean(inimoms[np.nonzero(inimoms)])
##            print(sta,moment,phs)
            if phs is 'P':
                outpfile.write('%.8e %s\n' % (moment,sta))
            elif phs is 'S':
                outsfile.write('%.8e %s\n' % (moment,sta))

        outpfile.close()
        outsfile.close()
