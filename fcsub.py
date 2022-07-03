import os
import subprocess
from antelope.datascope import *
import numpy as np
import math
import globaldb as g
from scipy.signal import *
import matplotlib.pyplot as plt
import pymutt
import seissign as seis
## Subrountines for inverting t*
## Written by S. Wei, June 2014
## readcoda: read S coda and calculate amplitude spectrum
## caldist: calculate the interval between two events


def readcoda(ordir,esinfo,sta,schan,winpara):
## Read S coda and calculate amplitude spectrum
## USAGE: (spec,freq,dospec,c_freq)=readcoda(
##                  ordir,esinfo,sta,schan,winpara)
    orid=ordir.split('_')[1]
    spec=np.zeros(100)
    freq=np.zeros(100)
    c_freq=np.zeros(100)
    dospec=False
    ## GET P TRAVEL TIME
    iphase='P'
##    g.dbja.record=g.dbja.find('iphase=~/%s/ && orid==%s && sta=~/%s/' % (iphase,orid,sta))
##    if g.dbja.record < 0:
##        print('No P arrival for %s %s'%(orid,sta))
##        return spec,freq,dospec,c_freq
    try:
        findresP=g.dbja.find('iphase=~/%s/ && orid==%s && sta=~/%s/' % (iphase,orid,sta))
    except DbfindEnd or DbfindError:
        print('No P arrival for %s %s'%(orid,sta))
        return spec,freq,dospec,c_freq
    else:
        g.dbja.record=findresP
    ptime=g.dbja.getv('time')[0]
    Ptt=ptime-esinfo['etime']
    ## GET S TRAVEL TIME
    iphase='S'
##    g.dbja.record=g.dbja.find('iphase=~/%s/ && orid==%s && sta=~/%s/' % (iphase,orid,sta))
##    if g.dbja.record < 0:
##        return spec,freq,dospec,c_freq
####        print('No S arrival for %s %s'%(orid,sta))
##        sss='''
##taup_time -ph s,S -h %f -evt %f %f -sta %f %f --time
##''' % (esinfo['edep'],esinfo['elat'],esinfo['elon'],esinfo['slat'],esinfo['slon'])
##        osout=os.popen(sss).read()
##        Stt=float(osout.split()[0])
####        print('TauP gives S travel time = %f' % Stt)
##    else:
##        stime=g.dbja.getv('time')[0]
##        Stt=stime-esinfo['etime']
    try:
        findresS=g.dbja.find('iphase=~/%s/ && orid==%s && sta=~/%s/' % (iphase,orid,sta))
    except DbfindEnd or DbfindError:
        return spec,freq,dospec,c_freq
##        print('No S arrival for %s %s'%(orid,sta))
        sss='''
taup_time -ph s,S -h %f -evt %f %f -sta %f %f --time
''' % (esinfo['edep'],esinfo['elat'],esinfo['elon'],esinfo['slat'],esinfo['slon'])
        osout=os.popen(sss).read()
        Stt=float(osout.split()[0])
##        print('TauP gives S travel time = %f' % Stt)
    else:
        g.dbja.record=findresS
        stime=g.dbja.getv('time')[0]
        Stt=stime-esinfo['etime']
       
    ## LOAD S CODA WAVE
    sacfl=ordir+'/'+orid+'_'+sta+'.'+schan+'.txt'
    if not os.path.isfile(sacfl):
        print('ERROR: %s does not exist' % sacfl)
        return spec,freq,dospec,c_freq
    dd=np.fromstring(("".join(open(sacfl).readlines()[30:])),sep=' ')
    dd=detrend(dd-np.mean(dd))
    dd=seis.seisfilter(dd,winpara['dt'],'highpass',winpara['lowf'])
    nn=dd.size
    tt=-winpara['pretime']+winpara['dt']*np.array(range(nn))
    tnbeg=-1-winpara['wln']
    tnend=-1
    tcbeg=winpara['startcoda']*Stt-Ptt
    tcend=tcbeg+winpara['wlc']
    nind=np.all((tt>=tnbeg,tt<=tnend),axis=0)
    noise=seis.taper(detrend(dd[nind]-np.mean(dd[nind])))
    cind=np.all((tt>=tcbeg,tt<=tcend),axis=0)
    coda=seis.taper(detrend(dd[cind]-np.mean(dd[cind])))
    ## CALCULATE SPECTURM
    (n_spec,n_freq)=seis.spec(noise,winpara['dt'],winpara['maxfreq'],'amp')
    (c_spec,c_freq)=seis.spec(coda,winpara['dt'],winpara['maxfreq'],'amp')    
    n_freq=n_freq[0:min(n_freq.size,c_freq.size)]
    c_freq=c_freq[0:min(n_freq.size,c_freq.size)]
    n_spec=n_spec[0:min(n_freq.size,c_freq.size)]
    c_spec=c_spec[0:min(n_freq.size,c_freq.size)]
    n_spec=seis.smooth(n_spec,window_len=winpara['smwin'])
    c_spec=seis.smooth(c_spec,window_len=winpara['smwin'])
    ## FIND THE LONGEST SEGMENTS OF CODA SPECTRA WITH SNR > SNRCRT
    snr=c_spec/n_spec
    lenspec=len(snr)
    snrcrt=winpara['snr'][0]
    w=0
    m=[]
    bindex=[]
    eindex=[]
    for kk in range(1,lenspec):
        if snr[kk]<snrcrt and snr[kk-1]>=snrcrt and kk==1:        # only first > crt
            w=1
            m.append(w)
            bindex.append(kk-w)
            eindex.append(kk-1)
            w=0
        elif snr[kk]>=snrcrt and snr[kk-1]>=snrcrt and kk==1:     # at first and continuously > crt
            w=w+2
        elif snr[kk]>=snrcrt and snr[kk-1]<snrcrt and kk>=1 and kk<(lenspec-1):    # begin of continuously > crt
            w=w+1
        elif snr[kk]>=snrcrt and snr[kk-1]>=snrcrt and kk>1 and kk<(lenspec-1):   # continuously >= crt
            w=w+1
        elif snr[kk]<snrcrt and snr[kk-1]>=snrcrt and kk>1 and kk<=(lenspec-1):    # end of continuously > crt
            m.append(w)
            bindex.append(kk-w)
            eindex.append(kk-1)
            w=0
        elif snr[kk]<snrcrt and snr[kk-1]<snrcrt and kk>=1 and kk<=(lenspec-1):     # continuously < crt
            w=0
        elif snr[kk]>=snrcrt and snr[kk]>=snrcrt and kk==(lenspec-1):     # at last and continuously > crt
            w=w+1
            m.append(w)
            bindex.append(kk-w+1)
            eindex.append(kk)
        elif snr[kk]>=snrcrt and snr[kk]<snrcrt and kk==(lenspec-1):      # only last > crt
            w=1
            m.append(w)
            bindex.append(kk-w+1)
            eindex.append(kk)
    if len(m)==0:
        return spec,freq,dospec,c_freq
    ## FIND THE LONGEST SEGMENT
    longest=m.index(max(m))
    frmin=c_freq[bindex[longest]]
    frmax=c_freq[eindex[longest]]
##    frangep=frmax-frmin
    frangep=np.log10(frmax)-np.log10(frmin)
    ## FAVOR THE SECOND LONGEST SEGMENT IF IT HAS LOWER FREQUENCY AND LONGER
    ## THAN 1/3 OF THE LONGEST ONE
    if len(m)>=2:
        longest2=m.index(sorted(m)[-2])
        frmin2=c_freq[bindex[longest2]]
        frmax2=c_freq[eindex[longest2]]
##        frangep2=frmax2-frmin2
        frangep2=np.log10(frmax2)-np.log10(frmin2)
        if frmin2<frmin and 3*frangep2>frangep and frangep2>winpara['snr'][1]:
            frmin=frmin2
            frmax=frmax2
            frangep=frangep2
            longest=longest2
    ## SPECTRUM OF P CACULATED FROM THE LONGEST OR THE SECOND LONGEST SEGMENT
    if frmax<10 and frangep<0.3:
        return spec,freq,dospec,c_freq
    elif frmax>10 and frangep<winpara['snr'][1]:
        return spec,freq,dospec,c_freq
    spec=c_spec[bindex[longest]:eindex[longest]]
    freq=c_freq[bindex[longest]:eindex[longest]]
    dospec=True
    
    ## PLOT WAVEFORM AND SPECTRA
    if winpara['plotspec']&dospec:
        plt.figure(1)
        plt.clf()
        plt.subplot(2,1,1)
        plt.xlabel(schan+' (seconds)')
        plt.ylabel('Velocity (nm/s)')
        plt.title('ORID: %s   Station: %s  S-coda: %.1f - %.1f s' % (ordir,sta,tcbeg,tcend))
        plt.plot(tt+Ptt,dd)
        plt.plot([tnbeg+Ptt,tnbeg+Ptt],[min(dd),max(dd)],'k',linewidth=2)
        plt.plot([tnend+Ptt,tnend+Ptt],[min(dd),max(dd)],'k',linewidth=2)
        plt.plot([tcend+Ptt,tcend+Ptt],[min(dd),max(dd)],'g',linewidth=2)
        plt.plot([tcbeg+Ptt,tcbeg+Ptt],[min(dd),max(dd)],'g',linewidth=2)
        plt.plot([Stt,Stt],[min(dd),max(dd)],'r',linewidth=2)
        plt.text(tnbeg+Ptt,min(dd),'Noise')
        plt.text(tcbeg+Ptt,min(dd),'S coda')
        plt.text(Stt,min(dd),'S arrival')
        plt.xlim([min(0,min(tt)+Ptt),min(max(tt)+Ptt,tcend+Ptt+20)])
        plt.subplot(2,1,2)
        plt.semilogy(c_freq,c_spec,'b')
        plt.semilogy(n_freq,n_spec,'r')
        plt.semilogy([frmin,frmin],[min(n_spec),max(c_spec)],'g')
        plt.semilogy([frmax,frmax],[min(n_spec),max(c_spec)],'g')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude Spectrum')
##        plt.show()
##        s=raw_input('Is the S coda good? (Y)')
##        if s.upper()=='N':
##            dospec=False
##        else:
##            dospec=True
        plt.savefig(winpara['figdir']+'%s_%s_Scoda.eps' % (ordir,sta))

    return spec,freq,dospec,c_freq



def caldist(lat1,lon1,dep1,lat2,lon2,dep2):
## Calculate the interval between two events
## USAGE: dist=caldist(lat1,lon1,dep1,lat2,lon2,dep2)
    Re=6371.0
    phi1=math.radians(90.0-lat1)
    theta1=math.radians(lon1)
    r1=Re-dep1
    x1=r1*math.cos(theta1)*math.sin(phi1)
    y1=r1*math.sin(theta1)*math.sin(phi1)
    z1=r1*math.cos(phi1)
    phi2=math.radians(90.0-lat2)
    theta2=math.radians(lon2)
    r2=Re-dep2
    x2=r2*math.cos(theta2)*math.sin(phi2)
    y2=r2*math.sin(theta2)*math.sin(phi2)
    z2=r2*math.cos(phi2)
    dist=math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    return dist


##def rmsenv(ordir,esinfo,sta,chan,winpara):
#### Read S coda and calculate RMS of narrow band envelopes
#### USAGE: fcsub.rmsenv(ordir,esinfo,sta,chan[1:3],winpara)
##    orid=ordir.split('_')[1]
##    spec=np.zeros(100)
##    freq=np.zeros(100)
##    c_freq=np.zeros(100)
##    dospec=False
##    ## GET P TRAVEL TIME
##    iphase='P'
##    g.dbja.record=g.dbja.find('iphase=~/%s/ && orid==%s && sta=~/%s/' % (iphase,orid,sta))
##    if g.dbja.record < 0:
##        print('No P arrival for %s %s'%(orid,sta))
##        return spec,freq,dospec,c_freq
##    ptime=g.dbja.getv('time')[0]
##    Ptt=ptime-esinfo['etime']
##    ## GET S TRAVEL TIME
##    iphase='S'
##    g.dbja.record=g.dbja.find('iphase=~/%s/ && orid==%s && sta=~/%s/' % (iphase,orid,sta))
##    if g.dbja.record < 0:
##        return spec,freq,dospec,c_freq
####        print('No S arrival for %s %s'%(orid,sta))
##        sss='''
##taup_time -ph s,S -h %f -evt %f %f -sta %f %f --time
##''' % (esinfo['edep'],esinfo['elat'],esinfo['elon'],esinfo['slat'],esinfo['slon'])
##        osout=os.popen(sss).read()
##        Stt=float(osout.split()[0])
####        print('TauP gives S travel time = %f' % Stt)
##    else:
##        stime=g.dbja.getv('time')[0]
##        Stt=stime-esinfo['etime']
##    ## LOAD EACH COMPONENT
##    for ichan in range(len(chan)):
##        sacfl=ordir+'/'+orid+'_'+sta+'.'+chan[ichan]+'.txt'
##        if not os.path.isfile(sacfl):
##            print('ERROR: %s does not exist' % sacfl)
##            return
##        dd=np.fromstring(("".join(open(sacfl).readlines()[30:])),sep=' ')
##        dd=detrend(dd-np.mean(dd))
##        for ifreq in range(len(winpara['freq'])):
##            dd=seis.seisfilter(dd,winpara['dt'],'bandpass',winpara['freq'][ifreq])
##            tnbeg=-1-winpara['wln']+Ptt
##            tnend=-1+Ptt
##            tcbeg=winpara['startcoda']*Stt
##            tcend=tcbeg+winpara['wlc']+Ptt
####            nind=np.all((tt>=tnbeg,tt<=tnend),axis=0)
####            noise=seis.taper(detrend(dd[nind]-np.mean(dd[nind])))
####            cind=np.all((tt>=tcbeg,tt<=tcend),axis=0)
####            coda=seis.taper(detrend(dd[cind]-np.mean(dd[cind])))
####            aind=np.all((tt>=0,tt<=200),axis=0)
####            if ichan==0:
####                wave=seis.taper(detrend(dd[aind]-np.mean(dd[aind])))
####            else:
####                wave=np.vstack((wave,seis.taper(detrend(dd[aind]-np.mean(dd[aind])))))
##            if ifreq==0:
##                waveblock=np.atleast_3d(dd)
##            else:
##                waveblock=np.dstack((waveblock,dd))
##        if ichan==0:
##            nn=dd.size
##            tt=-winpara['pretime']+winpara['dt']*np.array(range(nn))+Ptt
##            wave=waveblock
##        else:
##            wave=np.vstack((wave,waveblock))
##    rms=np.zeros((len(winpara['freq']),len(tt)))
##    temp=np.zeros(len(tt))
##    for ifreq in range(len(winpara['freq'])):
##        for itime in range(len(tt)):
##            temp[itime]=np.sqrt(np.mean((wave[:,itime,ifreq]**2)))
##        rms[ifreq]=seis.smooth(temp,int(10/winpara['dt']),window='flat')
##    plt.figure(2)
##    plt.clf()
##    for ifreq in range(len(winpara['freq'])):
##        plt.semilogy(tt,rms[ifreq])
##        plt.text(20,max(rms[ifreq]),'%.2f-%.2f (Hz)' % tuple(winpara['freq'][ifreq]))
##    plt.xlim([0,200])
####    plt.axvline(tnbeg,color='k',linewidth=2)
####    plt.axvline(tnend,color='k',linewidth=2)
##    plt.axvline(tcend,color='g',linewidth=2)
##    plt.axvline(tcbeg,color='g',linewidth=2)
####    plt.axvline(Stt,color='r',linewidth=2)
####    plt.text(tnbeg,rms.min(),'Noise')
##    plt.text(tcbeg,rms.min(),'S coda')
####    plt.text(Stt,rms.min(),'S arrival')
##    plt.savefig(winpara['figdir']+'%s_%s_RMS.eps' % (ordir,sta))
##    return


def readpspec(ordir,esinfo,sta,pchan,winpara):
## Read P and calculate amplitude spectrum
## USAGE: (spec,freq,dospec,c_freq)=readcoda(
##                  ordir,esinfo,sta,pchan,winpara)
    orid=ordir.split('_')[1]
    spec=np.zeros(100)
    freq=np.zeros(100)
    c_freq=np.zeros(100)
    dospec=False
    ## GET P TRAVEL TIME
    iphase='P'
##    g.dbja.record=g.dbja.find('iphase=~/%s/ && orid==%s && sta=~/%s/' % (iphase,orid,sta))
##    if g.dbja.record < 0:
##        print('No P arrival for %s %s'%(orid,sta))
##        return spec,freq,dospec,c_freq
    try:
        findresP=g.dbja.find('iphase=~/%s/ && orid==%s && sta=~/%s/' % (iphase,orid,sta))
    except DbfindEnd or DbfindError:
        print('No P arrival for %s %s' % (orid,sta))
        return spec,freq,dospec,c_freq
    else:
        g.dbja.record=findresP
    ptime=g.dbja.getv('time')[0]
    Ptt=ptime-esinfo['etime']
    ## LOAD P WAVE
    sacfl=ordir+'/'+orid+'_'+sta+'.'+pchan+'.txt'
    if not os.path.isfile(sacfl):
        print('ERROR: %s does not exist' % sacfl)
        return spec,freq,dospec,c_freq
    dd=np.fromstring(("".join(open(sacfl).readlines()[30:])),sep=' ')
    dd=detrend(dd-np.mean(dd))
    dd=seis.seisfilter(dd,winpara['dt'],'highpass',winpara['lowf'])
    nn=dd.size
    tt=-winpara['pretime']+winpara['dt']*np.array(range(nn))
    tnbeg=-0.5-winpara['wln']
    tnend=-0.5
    tcbeg=-0.5
    tcend=-0.5+winpara['wlc']
    nind=np.all((tt>=tnbeg,tt<=tnend),axis=0)
    noise=seis.taper(detrend(dd[nind]-np.mean(dd[nind])))
    cind=np.all((tt>=tcbeg,tt<=tcend),axis=0)
    pdata=seis.taper(detrend(dd[cind]-np.mean(dd[cind])))
    ## CALCULATE SPECTURM
    nft=1024
    npi=3.0
    mtmresult=pymutt.mtft(noise,dt=winpara['dt'],npi=npi,nwin=int(npi*2-1),paddedlen=nft)
    newspec=mtmresult['power'][1:]
    n_freq=(mtmresult['df']*np.arange(len(mtmresult['power'])))[1:]
    n_spec=np.sqrt(newspec)/(2*np.pi*n_freq)
    mtmresult=pymutt.mtft(pdata,dt=winpara['dt'],npi=npi,nwin=int(npi*2-1),paddedlen=nft)
    newspec=mtmresult['power'][1:]
    c_freq=(mtmresult['df']*np.arange(len(mtmresult['power'])))[1:]
    c_spec=np.sqrt(newspec)/(2*np.pi*n_freq)
    n_spec=seis.smooth(n_spec,window_len=winpara['smwin'])
    c_spec=seis.smooth(c_spec,window_len=winpara['smwin'])
    ## FIND THE LONGEST SEGMENTS OF CODA SPECTRA WITH SNR > SNRCRT
    snr=c_spec/n_spec
    lenspec=len([ifreq for ifreq in freq[0] if ifreq<15])
##    lenspec=len(snr)
    snrcrt=winpara['snr'][0]
    w=0
    m=[]
    bindex=[]
    eindex=[]
    for kk in range(1,lenspec):
        if snr[kk]<snrcrt and snr[kk-1]>=snrcrt and kk==1:        # only first > crt
            w=1
            m.append(w)
            bindex.append(kk-w)
            eindex.append(kk-1)
            w=0
        elif snr[kk]>=snrcrt and snr[kk-1]>=snrcrt and kk==1:     # at first and continuously > crt
            w=w+2
        elif snr[kk]>=snrcrt and snr[kk-1]<snrcrt and kk>=1 and kk<(lenspec-1):    # begin of continuously > crt
            w=w+1
        elif snr[kk]>=snrcrt and snr[kk-1]>=snrcrt and kk>1 and kk<(lenspec-1):   # continuously >= crt
            w=w+1
        elif snr[kk]<snrcrt and snr[kk-1]>=snrcrt and kk>1 and kk<=(lenspec-1):    # end of continuously > crt
            m.append(w)
            bindex.append(kk-w)
            eindex.append(kk-1)
            w=0
        elif snr[kk]<snrcrt and snr[kk-1]<snrcrt and kk>=1 and kk<=(lenspec-1):     # continuously < crt
            w=0
        elif snr[kk]>=snrcrt and snr[kk]>=snrcrt and kk==(lenspec-1):     # at last and continuously > crt
            w=w+1
            m.append(w)
            bindex.append(kk-w+1)
            eindex.append(kk)
        elif snr[kk]>=snrcrt and snr[kk]<snrcrt and kk==(lenspec-1):      # only last > crt
            w=1
            m.append(w)
            bindex.append(kk-w+1)
            eindex.append(kk)
    if len(m)==0:
        return spec,freq,dospec,c_freq
    ## FIND THE LONGEST SEGMENT
    longest=m.index(max(m))
    frmin=c_freq[bindex[longest]]
    frmax=c_freq[eindex[longest]]
##    frangep=frmax-frmin
    frangep=np.log10(frmax)-np.log10(frmin)
    ## FAVOR THE SECOND LONGEST SEGMENT IF IT HAS LOWER FREQUENCY AND LONGER
    ## THAN 1/3 OF THE LONGEST ONE
    if len(m)>=2:
        longest2=m.index(sorted(m)[-2])
        frmin2=c_freq[bindex[longest2]]
        frmax2=c_freq[eindex[longest2]]
##        frangep2=frmax2-frmin2
        frangep2=np.log10(frmax2)-np.log10(frmin2)
        if frmin2<frmin and 3*frangep2>frangep and frangep2>winpara['snr'][1]:
            frmin=frmin2
            frmax=frmax2
            frangep=frangep2
            longest=longest2
    ## SPECTRUM OF P CACULATED FROM THE LONGEST OR THE SECOND LONGEST SEGMENT
    if frmax<10 and frangep<0.3:
        return spec,freq,dospec,c_freq
    elif frmax>10 and frangep<winpara['snr'][1]:
        return spec,freq,dospec,c_freq
    spec=c_spec[bindex[longest]:eindex[longest]]
    freq=c_freq[bindex[longest]:eindex[longest]]
    c_freq=c_freq[0:lenspec]
    dospec=True
    
    ## PLOT WAVEFORM AND SPECTRA
    if winpara['plotspec']&dospec:
        plt.figure(1)
        plt.clf()
        plt.subplot(2,1,1)
        plt.xlabel(pchan+' (seconds)')
        plt.ylabel('Velocity (nm/s)')
        plt.title('ORID: %s   Station: %s  S-coda: %.1f - %.1f s' % (ordir,sta,tcbeg,tcend))
        plt.plot(tt+Ptt,dd)
        plt.plot([tnbeg+Ptt,tnbeg+Ptt],[min(dd),max(dd)],'k',linewidth=2)
        plt.plot([tnend+Ptt,tnend+Ptt],[min(dd),max(dd)],'k',linewidth=2)
        plt.plot([tcend+Ptt,tcend+Ptt],[min(dd),max(dd)],'g',linewidth=2)
        plt.plot([tcbeg+Ptt,tcbeg+Ptt],[min(dd),max(dd)],'g',linewidth=2)
        plt.text(tnbeg+Ptt,min(dd),'Noise')
        plt.text(tcbeg+Ptt,min(dd),'P')
        plt.xlim([tnbeg+Ptt-10,tcend+Ptt+10])
        plt.subplot(2,1,2)
        plt.semilogy(c_freq,c_spec,'b')
        plt.semilogy(n_freq,n_spec,'r')
        plt.semilogy([frmin,frmin],[min(n_spec),max(c_spec)],'g')
        plt.semilogy([frmax,frmax],[min(n_spec),max(c_spec)],'g')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude Spectrum')
        plt.savefig(winpara['figdir']+'%s_%s_Pspec.eps' % (ordir,sta))

    return spec,freq,dospec,c_freq
