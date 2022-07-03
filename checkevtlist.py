#!/opt/antelope/5.3/bin/python
## PROGRAM TO CHECK THE ERRORS OF THE EVENT LIST
## WRITTEN BY S. WEI, SEP. 2014

import os,sys
sys.path.append(os.environ['ANTELOPE'] +'/data/python')
from antelope.datascope import *
import globaldb as g

subdb='2'
## IMPORT GLOBAL VARIABLES
g.dbread(subdb)

## READ EVENT LIST
evtlst='/P/weisq/attentomo/input/checklst/sub'+subdb+'.lst'
outlst=open('/P/weisq/attentomo/input/eventid_sub'+subdb+'.lst','w')
orid1=[int(line.split()[0]) for line in open(evtlst).readlines()[1:]]
orid2=[int(line.split()[1]) for line in open(evtlst).readlines()[1:]]
mb=[line.split()[2] for line in open(evtlst).readlines()[1:]]
evid0=[int(line.split()[3]) for line in open(evtlst).readlines()[1:]]
lat0=[float(line.split()[4]) for line in open(evtlst).readlines()[1:]]
lon0=[float(line.split()[5]) for line in open(evtlst).readlines()[1:]]
dep0=[float(line.split()[6]) for line in open(evtlst).readlines()[1:]]

## CHECKING
for id in range(len(evid0)):
    g.dbass.record = g.dbass.find('orid==%d' % orid2[id])
    if g.dbass.record < 0:
        print('Error: Could not find arrival corresponding to orid2==%d ' % orid2[id])
        continue
    evla2 = g.dbass.getv('lat')[0]
    evlo2 = g.dbass.getv('lon')[0]
    evdp2 = g.dbass.getv('depth')[0]
    evid2 = g.dbass.getv('evid')[0]
    g.dbass.record = g.dbass.find('orid==%d' % orid1[id])
    if g.dbass.record < 0:
        print('Warning: Could not find arrival corresponding to orid1==%d ' % orid1[id])
        evid1=evid2
    else:
        evid1 = g.dbass.getv('evid')[0]
    if evid1==evid2 and evid0[id]==1111:
        pass
    elif evid1==evid2 and evid1==evid0[id]:
##    if evid1==evid2:
        pass
    else:
        print('Error: EVIDs do not match')
        print(orid1[id],orid2[id],evid0[id],lat0[id],lon0[id],dep0[id],evid1,evid2)
        continue
##    if abs(evla2-lat0[id])>0.1 or abs(evlo2-lon0[id])>0.1 or abs(evdp2-dep0[id])>10:
##        print('Error: locations do not match')
##        print(orid1[id],orid2[id],evid0[id],lat0[id],lon0[id],dep0[id],evla2,evlo2,evdp2)
    outlst.write('%d %s %.4f %.4f %.4f\n' %(orid2[id],mb[id],evla2,evlo2,evdp2))
