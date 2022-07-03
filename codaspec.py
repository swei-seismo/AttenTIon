#!/opt/antelope/5.4/bin/python
## LOAD S-CODA WAVE AND CALCULATE SPECTRUM
## WRITTEN BY S. WEI, JUNE 2014
## evt - ordir:%s - mb:%f
##                  time:%f
##                  lat:%f
##                  lon:%f
##                  dep:%f

import os,glob
import sys
sys.path.append(os.environ['ANTELOPE'] +'/data/python')
from antelope.datascope import *
import numpy as np
from scipy.fftpack import *
from scipy.signal import detrend
import matplotlib.pyplot as plt
import globaldb as g
import fcsub

sacdir='/P/weisq/attentomo/sacfl/'
fcdir='/P/weisq/attentomo/fc_1.5start30slogwin/'
specfigdir=fcdir+'specfig/'
specdir=fcdir+'spec/'
databases=[11,22,33,44]
if not os.path.isdir(fcdir):
    os.mkdir(fcdir)
if not os.path.isdir(specfigdir):
    os.mkdir(specfigdir)
if not os.path.isdir(specdir):
    os.mkdir(specdir)
orlstfl=open(specdir+'orid.lst','a')

## PARAMETERS FOR WINDOWING
winpara={}
winpara['pretime']=50
winpara['dt']=0.025
winpara['wln']=30
winpara['wlc']=30
winpara['plotspec']=False
winpara['figdir']=specfigdir
winpara['maxfreq']=15
winpara['snr']=[2.0,0.3]     ## [min SNR, min length of log10(frequency band)]
winpara['lowf']=0.02
winpara['smwin']=21
winpara['startcoda']=1.5
winpara['freq']=[[0.25,0.5],[0.5,1.0],[1.0,2.0],[2.0,4.0],[4.0,8.0],[8.0,16.0]]

## READ STATION INFOMATION
stafl='/P/weisq/attentomo/input/station.lst'
staloc={line.split()[0]:[float(line.split()[2]),float(line.split()[1])]
     for line in open(stafl).readlines()[1:]}

## PAIRING orid AND SUBDATABASE
evtall={}
for subdb in databases:
    evlst='/P/weisq/attentomo/input/eventid_sub%d.lst' % (subdb)
    oridlst=[line.split()[0] for line in open(evlst).readlines()[1:]]
    maglst=[float(line.split()[7]) for line in open(evlst).readlines()[1:]]
    for ior in range(len(oridlst)):
        orid=str(subdb)+'_'+oridlst[ior]
        evtall[orid]=maglst[ior]
os.chdir(sacdir)
evt={}
for ordir in glob.glob('*_*'):
    if evtall.keys().count(ordir)==1:
        evt[ordir]={}
        evt[ordir]['mb']=evtall[ordir]

evtcount1=0
evtcount2=0

## LOAD S CODA FROM ANTELOPE
for subdb in databases:
    print('Subdb = %d' % subdb)
    ## LOAD DATABASE AND GLOBAL VARIABLES
    g.dbread(str(subdb))
    ## READ EACH EVENT
    for ordir in evt:
##    for ordir in ['4_54417']:
        if int(ordir.split('_')[0])!=subdb:
            continue
        if os.path.isfile(specdir+'%s.stn' % ordir):
            print('%s.stn already exists' % ordir)
##            continue
        orid=ordir.split('_')[1]
        print('Read orid: %s' % ordir)
        evtcount1=evtcount1+1
        stnlstfl=open(specdir+'%s.stn' % ordir,'w')
        try:
            findres=g.dbor.find('orid==%s' % (orid))
        except DbfindEnd or DbfindError:
            print('Error: No orid %d' % orid)
            continue
        else:
            g.dbor.record=findres
##        g.dbor[3]=g.dbor.find('orid==%s' % (orid))
        (evt[ordir]['time'],evt[ordir]['lat'],evt[ordir]['lon'],evt[ordir]['dep']
             )=g.dbor.getv('time','lat','lon','depth')
        ## READ EACH STATION OF EACH EVENT
        stnlst=[]   ## STATIONS WITH USEFUL SPECTRA
        for sacfl in glob.glob(ordir+'/*HZ*.txt'):
            sta=sacfl.split('.')[0].split('_')[2]
            if sta=='F01' or sta=='F02W':
                continue
##            print('Read %s %s' %(orid,sta))
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
            ## READ S CODA AND CALCULATE ITS SPECTRUM
            esinfo={}
            esinfo['etime']=evt[ordir]['time']
            esinfo['elat']=evt[ordir]['lat']
            esinfo['elon']=evt[ordir]['lon']
            esinfo['edep']=evt[ordir]['dep']
            esinfo['slat']=staloc[sta][0]
            esinfo['slon']=staloc[sta][1]
            (spec,freq,dospec,c_freq)=fcsub.readcoda(ordir,esinfo,sta,chan[1],winpara)
##            fcsub.rmsenv(ordir,esinfo,sta,chan[1:3],winpara)
            ## OUTPUT CODA SPECTRUM
            if dospec:
                allfreqfl=specdir+'allfreq'
                if not os.path.isfile(allfreqfl):
                    np.savetxt(allfreqfl,c_freq,fmt='%.3f')
##                print('Useful coda from %f to %f Hz' % (freq[0],freq[freq.size-1]))
                stnlst.append(sta)
                stnlstfl.write(sta+'\n')
                specfl=specdir+'%s_%s.spec' % (ordir,sta)
                np.savetxt(specfl,np.vstack((freq,spec)).transpose(),fmt='%.3f')
        stnlstfl.close()
        if len(stnlst)>4:
            evtcount2=evtcount2+1
            orlstfl.write('%s  %.2f  %.2f  %.2f  %.1f\n' % 
                (ordir,evt[ordir]['lat'],evt[ordir]['lon'],evt[ordir]['dep'],
                evt[ordir]['mb']))
            print('%d/%d have good S-coda' % (evtcount2,evtcount1))
orlstfl.close()
