#!/opt/antelope/5.4/bin/python
## PROGRAM TO EXTRACT SAC FILES FOR t* INVERSION
## Written by S. Wei, MAR. 2014

import os,glob,shutil
import subprocess
import sys,time
sys.path.append(os.environ['ANTELOPE'] +'/data/python')
from antelope.datascope import *
import numpy as np
import globaldb as g


s=raw_input('Please input orid, subdb: ')
orid=int(s.split()[0])
subdb=s.split()[1]

##subdb='1'
##evlst='/P/weisq/attentomo/input/eventid_sub%s.lst' % subdb
####oridlst=[int(line.split()[0]) for line in open(evlst)]
##orid=51839

workdir='/P/weisq/attentomo/sacfl/'
orienfl='/P/weisq/attentomo/input/station.lst'
respdir='/P/sopac/mydatabases/lau-alldb/instrument-responses/'
if not os.path.isfile(orienfl):
    print('Cannot find station.lst')
    exit()
print(orid,subdb)
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

## IMPORT GLOBAL VARIABLES
print('============Working on ORID # %d of Subdb %s============' % (orid,subdb))
g.dbread(subdb)

## FREQUENCY DOMAIN PARAMETERS AND ALPHA
loco=0.05
hico=15
##TIME DOMAIN PARAMETERS
pretime=50
posttime=400
dt=0.025

# READING STATION ORIENTATIONS
phi={line.split()[0]:line.split()[5] for line in open(orienfl).readlines()[1:]}

ordir=workdir+subdb+'_'+str(orid)
##if os.path.isdir(ordir):
##    print('%d from %s exists and skipped' % (orid,subdb))
##    exit()
##    shutil.rmtree(ordir)
if not os.path.isdir(ordir):
    os.mkdir(ordir)

## DETERMINE NUMBER OF STATIONS AND RECORDS
## ALL P+S RECORDS LOCATED BY HUMAN
dbsubSTA=g.dbja.subset('orid==%d && auth!~/orbassoc_l/' % orid)
dbsubSTA=dbsubSTA.sort('sta')
numSTA=dbsubSTA.query('dbRECORD_COUNT')
if numSTA==0:
    print('Zero stations for ORID # %d' % orid)
    exit()
## ALL P ARRIVELS LOCATED BY HUMAN
dbsubP=g.dbja.subset('orid==%d && iphase=~/P/ && auth!~/orbassoc_l/' % orid)
dbsubP=dbsubP.sort('sta')
numP=dbsubP.query('dbRECORD_COUNT')
## ALL S ARRIVELS LOCATED BY HUMAN
dbsubS=g.dbja.subset('orid==%d && iphase=~/S/ && auth!~/orbassoc_l/' % orid)
dbsubS=dbsubS.sort('sta')
numS=dbsubS.query('dbRECORD_COUNT')

## GET EVENT INFOMATION
g.dbor.record=g.dbor.find('orid==%d' % (orid))
##year=str(g.dbor.getv('jdate')[0])[0:4]
##jday=str(g.dbor.getv('jdate')[0])[4:7]

###################### LOOP OVER RECORDS IN SUBSET TABLE ####################
dbsub=dbsubP
for ii in range(numP):
##for ii in [47]:
    dbsub.record=ii
    sta=dbsub.getv('sta')[0]
    if sta=='AFI':
        print('Skip AFI')
        continue
    if sta!='MSVF':
##        print('Skip MSVF')
        continue

    ## DETERMINE THE CHANNEL NAMES
    if sta=='AFI':
        chan=['BHE_00','BHN_00','BHZ_00']
    elif sta=='MSVF':
        chan=['BH1_00','BH2_00','BHZ_00']
    elif len(sta)==4 and sta.endswith('W'): ## WHOI OBS
        chan=['HH1','HH2','HHZ']
    elif len(sta)==3 and not (sta == 'FOA'):  ## LDEO OBS
        chan=['BH2','BH1','BHZ']
    else:                                   ## Other land station
        chan=['BHE','BHN','BHZ']

##    if os.path.isfile(ordir+'/'+str(orid)+'_'+sta+'.'+chan[2]+'.txt'):
####        print('%d at %s exists, skip.' % (orid,sta))
##        continue
    print('Working on station %s  (%d of %d)' % (sta,ii+1,numP))

    if phi[sta]=='NaN':
        rotation=False
        phiangle=0
    else:
        rotation=True
        phiangle=float(phi[sta])
            
    iphase='P'

    ## GET EVENT ARRIVAL TIME, EVENT INFORMATION, AND STATION INFO
    try:
        findres=g.dbja.find('iphase=~/%s/ && orid==%d && sta=~/%s/' % (iphase,orid,sta))
    except DbfindEnd or DbfindError:
        print('Could not find arrival corresponding to iphase=~/%s/ && orid==%d && sta=~/%s/' % (iphase,orid,sta))
        continue
    else:
        g.dbja.record=findres
##    g.dbja.record=g.dbja.find('iphase=~/%s/ && orid==%d && sta=~/%s/' % (iphase,orid,sta))
##    if g.dbja.record < 0:
##        print('Could not find arrival corresponding to iphase=~/%s/ && orid==%d && sta=~/%s/' % (iphase,orid,sta))
    at  = g.dbja.getv('time')[0]
    ori = g.dbja.getv('orid')[0]
    g.dbass.record = g.dbass.find('orid==%d' % ori)
    evla = g.dbass.getv('lat')[0]
    evlo = g.dbass.getv('lon')[0]
    evdp = g.dbass.getv('depth')[0]
    baz  = g.dbass.getv('seaz')[0]
    g.dbresp.record=g.dbresp.find('sta=~/%s/' % sta)
    stla = g.dbresp.getv('lat')[0]
    stlo = g.dbresp.getv('lon')[0]
    atime = time.strftime('%Y %j %H %M %S', time.gmtime(at))
    [year,jday,hour,minute,sec]=atime.split()

    ## EXTRACTING SEISMOGRAM dd AND tt
    t0=at-pretime
    t1=at+posttime
    goodchan=[True,True,True]
    for ichan in range(len(chan)):
        tr=g.dbwf.trloadchan(t0,t1,sta,chan[ichan])
##        tr=trloadchan(g.dbwf,t0,t1,sta,chan[ichan])
        tr.record=0
        rawdata=tr.trdata()
        ## CHECK IF THE CHANNEL HAS DATA
        if np.mean(np.array(rawdata))==0 and np.median(np.array(rawdata))==0:
            goodchan[ichan]=False
        if ichan==0:
            rawdt=1.0/(tr.getv('samprate')[0])
            nn=int((pretime+posttime)/rawdt)
            if nn>=len(rawdata):
                rawdata=rawdata+(0,)*(nn-len(rawdata))
            else:
                rawdata=rawdata[0:nn]
            rawdd=np.array([rawdata])
            rawtt=(t0-at)+rawdt*np.array(range(nn))
        else:
            if nn>=len(rawdata):
                rawdata=rawdata+(0,)*(nn-len(rawdata))
            else:
                rawdata=rawdata[0:nn]
            rawdd=np.vstack((rawdd,rawdata))
        np.savetxt('%s.asc' % chan[ichan],
                   np.vstack((rawtt,rawdd[ichan])).transpose(), fmt='%1.7e')
        ## REMOVE INSTRUMENT RESPONSE
        if sta=='MSVF':
            respfl=respdir+'RESP.II.MSVF.00.BHZ'
        else:
            respfl=glob.glob(respdir+'RESP*'+sta+'*'+chan[ichan])[0]
        netid=(os.path.basename(respfl)).split('.')[1]
        locid=(os.path.basename(respfl)).split('.')[3]
        if locid=='':
            locid='*'
##        cmd = '''evalresp %s %s %s %s 0.001 %f 2049 -u 'vel' -f %s''' % (
##            sta,chan[ichan][0:3],year,jday,0.5/dt+5,respfl)
##        os.system(cmd)
##        afile=glob.glob('AMP*')[0]
##        pfile=glob.glob('PHASE*')[0]
        if ichan==0:    ## HH1 FOR WHOI, BH2 FOR LEDO, BHE FOR LAND
            cmpinc=90
            cmpaz=90-phiangle
        elif ichan==1:  ## HH2 FOR WHOI, BH1 FOR LEDO, BHN FOR LAND
            cmpinc=90
            cmpaz=0-phiangle
        elif ichan==2:  ## VERTICAL
            if len(sta)==3 and not (sta == 'FOA'):
                cmpinc=180
            else:
                cmpinc=0
            cmpaz=0
        if cmpaz<0:
            cmpaz=360+cmpaz
        print(sta,locid)
        sss = '''
readtable content xy %(ch)s.asc
w %(newfn)s.sac
r %(newfn)s.sac
ch stla %(sla)f
ch stlo %(slo)f
ch evla %(ela)f
ch evlo %(elo)f
ch evdp %(edp)f
ch cmpinc %(inc)d
ch cmpaz %(az)d
ch lovrok true
ch lcalda true
ch knetwk %(net)s
ch khole %(locid)s
ch kcmpnm %(ch2)s
ch nzyear %(year)s
ch nzjday %(jday)s
ch nzhour %(hour)s
ch nzmin %(min)s
ch nzsec %(sec)d
ch nzmsec %(msec)d
ch iztype ib
ch kevnm %(evt)s
interpolate delta %(dtt)f
w over
r %(newfn)s.sac
rtr
rmean
taper
transfer from evalresp fname %(rfl)s to vel freqlimits 0.005 0.01 %(f3)f %(f4)f
w over
quit
        ''' %({'newfn':(ordir+'/'+str(orid)+'_'+sta+'.'+chan[ichan]),
               'inc':cmpinc,'az':cmpaz,'f3':(0.25/dt),'f4':(0.50/dt),
               'dtt':dt,'rfl':respfl,'evt':subdb+str(orid),
               'sla':stla,'slo':stlo,'ela':evla,'elo':evlo,'edp':evdp,
               'ch':chan[ichan],'ch2':chan[ichan][0:3],
               'year':year,'jday':jday,'hour':hour,'min':minute,
               'sec':int(divmod(float(sec),1)[0],),
               'msec':int(divmod(float(sec),1)[1]*1000),
               'net':netid,'locid':locid,})
        p=subprocess.Popen(['sac'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT )
        sacout=p.communicate(sss)
        if sacout[0].find('Using')>0:
            print(sacout[0])
    os.system('rm *.asc')
    if not goodchan[2]:
        continue
    ## SAVE THE NONE ZERO HORIZONTAL CHANNEL AS NORTH-SOUTH CHANNEL
    if goodchan[0] and goodchan[1]:
        pass
    elif goodchan[1]:
        rotation=False
    elif goodchan[0]:
        rotation=False
        shutil.copyfile(ordir+'/'+str(orid)+'_'+sta+'.'+chan[0],
                        ordir+'/'+str(orid)+'_'+sta+'.'+chan[1])
    else:
        continue
    ## ROTATE TO RADIAL AND TRANSVERSE IF ORIENTATION KNOWN
    ## AND SAVE TRANSVERSE COMPONENT AS NORTH-SOUTH CHANNEL
    if rotation:
        sss = '''
r %(echan)s.sac %(nchan)s.sac
rotate to gcp
w over
q
        ''' % {'echan':(ordir+'/'+str(orid)+'_'+sta+'.'+chan[0]),
               'nchan':(ordir+'/'+str(orid)+'_'+sta+'.'+chan[1])}
        p=subprocess.Popen(['sac'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT )
        p.communicate(sss)
    ## FILTERING AND READ SAC FILES
    sss = '''
r %(echan)s.sac %(nchan)s.sac %(zchan)s.sac
w alpha %(echan)s.txt %(nchan)s.txt %(zchan)s.txt
q
    ''' % ({'echan':(ordir+'/'+str(orid)+'_'+sta+'.'+chan[0]),
           'nchan':(ordir+'/'+str(orid)+'_'+sta+'.'+chan[1]),
           'zchan':(ordir+'/'+str(orid)+'_'+sta+'.'+chan[2]),
           'lco':loco,'hco':hico})
    p=subprocess.Popen(['sac'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT )
    p.communicate(sss)

##os.system('rm AMP* PHASE*')
