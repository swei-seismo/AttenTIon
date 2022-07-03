import os
import subprocess
from antelope.datascope import *
import numpy as np
import globaldb as g
from scipy.signal import *
from scipy.linalg import lstsq
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.size']=16
mpl.rcParams['axes.formatter.limits']=[-2,2]
import pymutt
import seissign as seis
## Subrountines for inverting t*
## Written by S. Wei, Nov. 2013
## Edited by S. Wei for t*(S), Feb. 2014
## readseismo:  extract raw data from antelope database
## getdbinfo:   get informatio of station and event
## fixwin:      window seismic data with a fixed length
## dospec:      calculate spectra of whole data, P/S signal, and P/S noise

def readseismo(pretime,dt,subdb,orid,sta,chan):
## Read seismogram from extsac_tstar.py 
## USAGE: (dd,time,flag) = readseismo(pretime,dt,orid,sta,chan)
## INPUT:  (pretime)  ==> Seconds before arrival to subset
##         (dt)       ==> Scalar: 1/(samplerate of data)
##         (subdb)    ==> Subset of database
##         (orid)     ==> Origin ID of event
##         (sta)      ==> Station name
##         (chan)     ==> Channel names (i.e. ['BHE','BHN','BHZ'])
## OUTPUT: (dd)       ==> (3Xn) numpy.ndarray: [(E_CHAN),(N_CHAN),(Z_CHAN)]
##         (time)     ==> (1Xn) numpy.ndarray: (Relative Time)
##         (flag)     ==> If reading successful
    for ichan in range(len(chan)):
        sacfl=g.sacdir+subdb+'_'+str(orid)+'/'+str(orid)+'_'+sta+'.'+chan[ichan]+'.txt'
##        print(sacfl)
        if not os.path.isfile(sacfl):
            print('ERROR: %s does not exist' % sacfl)
            dd=np.array(range(18000))
            tt=np.array(range(18000))
            flag=False
            return dd,tt,flag
        ddchan=np.fromstring(("".join(open(sacfl).readlines()[30:])),sep=' ')
##        ddchan=ddchan/1e9
##        ddchan=np.genfromtxt(sacfl,skiprows=30,invalid_raise=False)
        nn=ddchan.size
##        ddchan=ddchan.reshape(1,nn)
        if ichan==0:
            dd=ddchan
        else:
            dd=np.vstack((dd,ddchan))
    flag=True
    tt=-pretime+dt*np.array(range(nn))
    
    return dd,tt,flag


def getdbinfo(sta,orid,ORIG):
## Get desired information about given database
## USAGE: (ARRIV,Ptt,Stt) = getdbinfo(sta,orid,ORIG)
## INPUT: (sta)  = station name
##        (orid) = origin id (orid)
##        (ORIG) = Origin Disctionary with Keys:
##            [orid] = origin id (orid)
##            [time] = origin time (unix seconds)
##            [lat]  = hypocenter latitude
##            [lon]  = hypocenter longitude
##            [dep]  = hypocenter depth
##            [ml]   = local magnitude
## OUTPUT: (ARRIV) = Arrival Dictionary with Keys:
##            [sta]   = station name
##            [ptime] = P arrival time (unix seconds)
##            [parid] = P arrival id (arid)
##            [pchan] = channel that P was picked on
##            [stime] = S arrival time (unix seconds)
##            [sarid] = S arrival id (arid)
##            [schan] = channel that S was picked on
##            [delta] = distance in degrees from sta to origin
##            [baz]   = station to event azimuth
##         (Ptt) = P travel time
##         (Stt) = S travel time
##       (SDATA) = If S arrival exists
    
    ARRIV={}
    ## GET P ARRIVEL INFORMATION
    ARRIV['sta']=sta
    iphase='P'
    try:
        findresP=g.dbja.find('iphase=~/%s/ && orid==%d && sta=~/%s/' % (iphase,orid,sta))
    except DbfindEnd or DbfindError:
        print('Could not find assoc-arrival corresponding to iphase=~/%s/ && orid==%d && sta=~/%s/' % (iphase,orid,sta))
        return
    else:
        g.dbja.record=findresP
    (ARRIV['ptime'],ARRIV['parid'],ARRIV['pchan'],ARRIV['delta'],ARRIV['baz']
         )=g.dbja.getv('time','arid','chan','delta','seaz')
    Ptt=ARRIV['ptime']-ORIG['time']
    ## GET S ARRIVEL INFORMATION
    iphase='S'
    try:
        findresS=g.dbja.find('iphase=~/%s/ && orid==%d && sta=~/%s/' % (iphase,orid,sta))
    except DbfindEnd or DbfindError:
        ARRIV['stime']=111
        ARRIV['sarid']=111
        if sta=='AFI' or sta=='MSVF':
            ARRIV['schan']='EEE_00'
        else:
            ARRIV['schan']='EEE'
        Stt=Ptt+1
        SDATA=False
    else:
        g.dbja.record=findresS
        (ARRIV['stime'],ARRIV['sarid'],ARRIV['schan'])=g.dbja.getv('time','arid','chan')
        Stt=ARRIV['stime']-ORIG['time']
        SDATA=True
    return ARRIV,Ptt,Stt,SDATA



def fixwin(dd,tt,dt,chan,ARRIV,prewin,WLP,WLS,SDATA,doplot,orid,sta):
## Window seismic data
## USAGE: (p_dd,s_dd,pn_dd,sn_dd,NOS2)
##            =tstarsub.fixwin(dd,tt,dt,ARRIV[ii+1],prewin,WLP,WLS,NOS,doplot,orid,sta)
## INPUT: (dd)      ==> (3Xn) numpy.ndarray: [(E_CHAN),(N_CHAN),(Z_CHAN)]
##        (tt)      ==> (1Xn) numpy.ndarray: (Relative Time)
##        (dt)      ==> Scalar: 1/(samplerate of data)
##        (chan)    ==> List of channal names
##        (ARRIV)   ==> Arrival Dictionary with Keys (see getdbinfo)
##        (prewin)  ==> Seconds before P[0]/S[1] arrival for windowing
##        (WLP)     ==> Window lenght for P in seconds
##        (WLS)     ==> Window lenght for S in seconds
##        (SDATA)   ==> Existence of S arrival
##        (doplot)  ==> Bool variable for plotting spectrum
##        (orid)    ==> origin id
##        (sta)     ==> station name
## OUTPUT: (p_dd)   ==> Windowed P data, starting prewin seconds before P arrival
##         (s_dd)   ==> Windowed S data
##         (pn_dd)  ==> P noise data, WL length ending prewin seconds before P arrival
##         (sn_dd)  ==> S noise data
##         (pd_dd)  ==> P coda data, right before S arrivel

    gaps=0.5      ## seconds between S noise and signal
    ## WINDOWING P
    pch=chan.index(ARRIV['pchan'])
    pind=np.all([(tt>=-1*prewin[0]),(tt<=(WLP-prewin[0]))],axis=0)
    p_dd=detrend(dd[pch][pind]-np.mean(dd[pch][pind]))
    pnind=np.all([(tt<=-1*prewin[0]),(tt>=(-WLP-prewin[0]))],axis=0)
    pn_dd=detrend(dd[pch][pnind]-np.mean(dd[pch][pnind]))
    ## MAKING SURE P WAVE AND NOISE HAVE SAME SIZE
    if p_dd.size>pn_dd.size:
        p_dd=p_dd[0:(pn_dd.size)]
    elif p_dd.size<pn_dd.size:
        pn_dd=s_dd[0:(p_dd.size)]
    if doplot:
        scalefac=1e3
        plt.figure(3)
        plt.clf()
        plt.subplot(2,1,2)
        plt.xlabel(ARRIV['pchan']+' (seconds)')
        plt.ylabel('Velocity Amplitude (nm/s)')
        plt.title('Station: %s' % sta)
#        plt.plot(tt,dd[pch])
        tmpdd=dd[pch]
        b, a = butter(4, 0.01, btype='highpass')
        filtmpdd = filtfilt(b, a, tmpdd)
        plt.plot(tt,filtmpdd)
#        plt.xlim([-WLP-prewin[0]-3,WLP-prewin[0]+3])
        plt.xlim([-10,10])
        plt.ylim([-7e4,7e4]);
#        plt.ylim([np.floor(min(filtmpdd)/scalefac)*scalefac,np.ceil(max(filtmpdd)/scalefac)*scalefac])
        plt.axvline(-1*prewin[0],color='g')
        plt.axvline(WLP-prewin[0],color='g')
        plt.axvline(-WLP-prewin[0],color='g')
        plt.text(-WLP-prewin[0],min(p_dd),'Noise')
        plt.text(-1*prewin[0],min(p_dd),'P wave')
        plt.savefig(g.figdir+'/%d_%s_Pwindata.pdf' % (orid,sta))
    ## WINDOWING S ON BOTH HORIZONTAL CHANNELS
    if SDATA:
        sminusp=ARRIV['stime']-ARRIV['ptime']
        if sminusp<(WLS+WLP+prewin[0]-prewin[1]+gaps):
            print('P & S arrivels are too close - proceeding as if no S pick')
            SDATA=False
    if SDATA:
##        sch=chan.index(ARRIV['schan'])
        ## HORIZONTAL CHANNEL 1
        sch=0
        sind=np.all([(tt>=(sminusp-prewin[1])),(tt<=(sminusp+WLS-prewin[1]))],axis=0)
        s_dd1=detrend(dd[sch][sind]-np.mean(dd[sch][sind]))
######## Noise defined as gaps s before S arrival
        snind=np.all([(tt<=(sminusp-prewin[1]-gaps)),(tt>=(sminusp-WLS-prewin[1]-gaps))],axis=0)
        sn_dd1=detrend(dd[sch][snind]-np.mean(dd[sch][snind]))
########## Noise defined as right before P arrival
##        snind=np.all([(tt<=-1*prewin[1]),(tt>=(-WLS-prewin[1]))],axis=0)
##        sn_dd=detrend(dd[sch][snind]-np.mean(dd[sch][snind]))
######## P coda defined as right before S arrival
        pcind=np.all([(tt<=(sminusp-prewin[1])),(tt>=(sminusp-WLS-prewin[1]))],axis=0)
        pc_dd1=detrend(dd[sch][pcind]-np.mean(dd[sch][pcind]))
        ## MAKING SURE S WAVE, NOISE AND P CODA HAVE SAME SIZE
        minlen=min(s_dd1.size,sn_dd1.size,pc_dd1.size)
        s_dd1=s_dd1[0:minlen]
        sn_dd1=sn_dd1[0:minlen]
        pc_dd1=pc_dd1[0:minlen]
        ## HORIZONTAL CHANNEL 2
        sch=1
        sind=np.all([(tt>=(sminusp-prewin[1])),(tt<=(sminusp+WLS-prewin[1]))],axis=0)
        s_dd2=detrend(dd[sch][sind]-np.mean(dd[sch][sind]))
######## Noise defined as right before S arrival
        snind=np.all([(tt<=(sminusp-prewin[1]-gaps)),(tt>=(sminusp-WLS-prewin[1]-gaps))],axis=0)
        sn_dd2=detrend(dd[sch][snind]-np.mean(dd[sch][snind]))
########## Noise defined as right before P arrival
##        snind=np.all([(tt<=-1*prewin[1]),(tt>=(-WLS-prewin[1]))],axis=0)
##        sn_dd=detrend(dd[sch][snind]-np.mean(dd[sch][snind]))
######## P coda defined as right before S arrival
        pcind=np.all([(tt<=(sminusp-prewin[1])),(tt>=(sminusp-WLS-prewin[1]))],axis=0)
        pc_dd2=detrend(dd[sch][pcind]-np.mean(dd[sch][pcind]))
        ## MAKING SURE S WAVE, NOISE AND P CODA HAVE SAME SIZE
        minlen=min(s_dd2.size,sn_dd2.size,pc_dd2.size)
        s_dd2=s_dd2[0:minlen]
        sn_dd2=sn_dd2[0:minlen]
        pc_dd2=pc_dd2[0:minlen]
        if doplot:
            scalefac=2e3
            plt.figure(4)
            plt.clf()
            plt.subplot(2,1,1)
#            plt.xlabel(ARRIV['schan']+' (seconds)')
            plt.xlabel(chan[0]+' (seconds)')
            plt.ylabel('Velocity Amplitude (nm/s)')
            plt.title('Station: %s' % sta)
#            plt.plot(tt,dd[0])
            tmpdd=dd[0]
            b, a = butter(4, 0.01, btype='highpass')
            filtmpdd = filtfilt(b, a, tmpdd)
            plt.plot(tt,filtmpdd)
########## P coda defined as right before S arrival
##            plt.axvline(sminusp-WLS-prewin[1],color='g')
##            plt.text(sminusp-WLP-prewin[1],min(s_dd),'P coda')
######## Noise defined as right before S arrival
##            plt.xlim([sminusp-WLS-prewin[1]-10,sminusp+WLS-prewin[1]+10])
##            plt.plot([sminusp-prewin[1],sminusp-prewin[1]],[min(dd[sch]),max(dd[sch])],'g')
##            plt.plot([sminusp+WLS-prewin[1],sminusp+WLS-prewin[1]],[min(dd[sch]),max(dd[sch])],'g')
##            plt.plot([sminusp-WLS-prewin[1],sminusp-WLS-prewin[1]],[min(dd[sch]),max(dd[sch])],'g')
#            plt.xlim([-WLP-prewin[1]-5,sminusp+WLS-prewin[1]+5])
#            plt.xlim([sminusp-WLS-gaps-prewin[1]-3,sminusp+WLS-prewin[1]+3])
            plt.xlim([sminusp-10,sminusp+10])
            plt.ylim([-5e5,6e5]);
#            plt.ylim([np.floor(min(filtmpdd)/scalefac)*scalefac,np.ceil(max(filtmpdd)/scalefac)*scalefac])
            plt.axvline(sminusp-prewin[1],color='g')
            plt.axvline(sminusp+WLS-prewin[1],color='g')
            plt.axvline(sminusp-prewin[1]-gaps,color='g')
            plt.axvline(sminusp-WLS-prewin[1]-gaps,color='g')
            plt.axvline(-1*prewin[1],color='g')
            plt.text(-prewin[1],min(s_dd1),'P arrival')
            plt.text(sminusp-WLS-prewin[1]-gaps,min(s_dd1),'Noise')
            plt.text(sminusp-prewin[1],min(s_dd1),'S wave')
########## Noise defined as right before P arrival
##            plt.xlim([-WLP-prewin[1]-10,sminusp+WLS-prewin[1]+10])
##            plt.axvline(sminusp-prewin[1],color='g')
##            plt.axvline(sminusp+WLS-prewin[1],color='g')
##            plt.axvline(-1*prewin[1],color='g')
##            plt.axvline(-WLS-prewin[1],color='g')
##            plt.text(-WLP-prewin[1],min(s_dd),'Noise')
##            plt.text(-prewin[1],min(s_dd),'P arrival')
##            plt.text(sminusp-prewin[1],min(s_dd),'S wave')
            plt.subplot(2,1,2)
#            plt.xlabel(ARRIV['schan']+' (seconds)')
            plt.xlabel(chan[1]+' (seconds)')
            plt.ylabel('Velocity Amplitude (nm/s)')
#            plt.plot(tt,dd[1])
            tmpdd=dd[1]
            b, a = butter(4, 0.01, btype='highpass')
            filtmpdd = filtfilt(b, a, tmpdd)
            plt.plot(tt,filtmpdd)
##            plt.xlim([-WLP-prewin[1]-10,sminusp+WLS-prewin[1]+10])
#            plt.xlim([sminusp-WLS-gaps-prewin[1]-3,sminusp+WLS-prewin[1]+3])
            plt.xlim([sminusp-10,sminusp+10])
            plt.ylim([-5e5,6e5]);
#            plt.ylim([np.floor(min(filtmpdd)/scalefac)*scalefac,np.ceil(max(filtmpdd)/scalefac)*scalefac])
            plt.axvline(sminusp-prewin[1],color='g')
            plt.axvline(sminusp+WLS-prewin[1],color='g')
            plt.axvline(sminusp-prewin[1]-gaps,color='g')
            plt.axvline(sminusp-WLS-prewin[1]-gaps,color='g')
#            plt.axvline(-1*prewin[1],color='g')
#            plt.text(-prewin[1],min(s_dd2),'P arrival')
            plt.text(sminusp-WLS-prewin[1]-gaps,min(s_dd2),'Noise')
            plt.text(sminusp-prewin[1],min(s_dd2),'S wave')
            plt.savefig(g.figdir+'/%d_%s_Swindata.pdf' % (orid,sta))
    else:
        s_dd1=p_dd
        sn_dd1=pn_dd
        pc_dd1=p_dd
        s_dd2=p_dd
        sn_dd2=pn_dd
        pc_dd2=p_dd
    return p_dd,pn_dd,s_dd1,sn_dd1,pc_dd1,s_dd2,sn_dd2,pc_dd2,SDATA


def longseg(snr,snrcrtpara,freq,minf=0.05,maxf=15.0):
## FIND THE LONGEST SEGMENT OF SPECTRA WITH SNR > SNRCRT
## USAGE: (begind,endind,frmin,frmax,frange)=longseg(snr,snrcrtp,freq)
## INPUT: snr        = SIGNAL-TO-NOISE RATIO
##        snrcrtpara = [MINIMUM SIGNAL-TO-NOISE, MINIMUM LENGTH OF SEGMENT]
##        freq       = FREQUENCY
## OUTPUT: begind = INDEX OF THE BEGIN OF THE LONGEST SEGMENT
##         endind = INDEX OF THE END OF THE LONGEST SEGMENT
##         frmin  = MINIMUM FREQUENCY OF THE LONGEST SEGMENT
##         frmax  = MAXIMUM FREQUENCY OF THE LONGEST SEGMENT
##         frange = frmax - frmin
    ## TAKE SPECTRUM < maxf (15 Hz) and > minf (0.1 Hz)
##    print('Find longest frequency band')
    lenspec=len([ifreq for ifreq in freq if (ifreq<maxf and ifreq>=minf)])
    ind1=int(min(np.nonzero(freq>=minf)[0]))
    ind2=int(max(np.nonzero(freq<maxf)[0]))
    w=0
    m=[]
    bindex=[]
    eindex=[]
    snrcrt=snrcrtpara[0]
    for kk in range(ind1+1,lenspec):
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
        frange=0
        frmin=6
        frmax=6
        begind=0
        endind=0
        return begind,endind,frmin,frmax,frange
    ## FIND THE LONGEST SEGMENT
    longest=m.index(max(m))
    frmin=freq[bindex[longest]]
    frmax=freq[eindex[longest]]
    frange=frmax-frmin
##    print(frmin,frmax,frange,m)
    ## FAVOR THE SECOND LONGEST SEGMENT IF IT HAS LOWER FREQUENCY < 4 Hz AND LONGER
    ## THAN 1/4 OF THE LONGEST ONE
    if len(m)>=2:
        for mind in list(reversed(range(len(m)))):
            mii=mind-len(m)
            longest2=m.index(sorted(m)[mii])
            frmin2=freq[bindex[longest2]]
            if frmin2<=2.0:
                frmax2=freq[eindex[longest2]]
                frange2=frmax2-frmin2
##                print(frmin2,frmax2,frange2,snrcrtpara[1])
                if frmin2<frmin and 4*frange2>frange and frange2>snrcrtpara[1]:
                    frmin=frmin2
                    frmax=frmax2
                    frange=frange2
                    longest=longest2
                    break
    begind=bindex[longest]
    endind=eindex[longest]
    ## EXTEND FREQUENCY BAND TO lowSNR
    if snrcrtpara[2]<snrcrtpara[0]:
        if begind>ind1+1:
            while snr[begind-1]<snr[begind] and snr[begind-1]>snrcrtpara[2] and begind-1>ind1+1:
                begind=begind-1
        if endind<ind2-1:
            while snr[endind+1]<snr[endind] and snr[endind+1]>snrcrtpara[2] and endind+1<ind2-1:
                endind=endind+1
        frmin=freq[begind]
        frmax=freq[endind]
        frange=frmax-frmin
    return begind,endind,frmin,frmax,frange


def dospec(pwindata,swindata1,swindata2,dt,SDATA,orid,sta,snrcrtp,snrcrts,
           linresid,chan,doplot):
## Calculate amplitude spectrum of windowed waveform using multi-taper method
## USAGE: (spec_px,freq_px,spec_sx,freq_sx,spec,freq,n_spec,n_freq,frmn,frmx,
##         goodP1,goodS1)=tstarsub.dospec(PWINDATA,SWINDATA1,SWINDATA2,dt,
##                      SDATA,orid,sta,snrcrtp,snrcrts,lincor,chan,doplot)
## INPUT:   pwindata  = P windowed data [0] and noise [1]
##          swindata1 = S windowed data [0] and noise [1] and P coda [2] on channel 1
##          swindata2 = S windowed data [0] and noise [1] and P coda [2] on channel 2
##          dt        = 1/sample rate
##          SDATA     = existence of S arrival (Bool variable)
##          orid      = origin id
##          sta       = station name
##          snrcrtp   = minimum [SNR,width,lowSNR] of good P in freqency domain
##          snrcrts   = minimum [SNR,width,lowSNR] of good S in freqency domain
##          lincor    = MINIMUM LINEAR CORRELATION COEFFICIENTS
##          chan      = CHANNELS
##          doplot    = Bool variable for plotting spectrum
## OUTPUT:  spec_px = spectrum of good P signal
##          freq_px = freqency range of spectrum of good P signal
##          spec_sx = spectrum of good S signal
##          freq_sx = freqency range of spectrum of good S signal
##          spec    = spectrum of all signal
##          freq    = freqency range of spectrum of all signal
##          n_spec  = fspectrum of all noise
##          n_freq  = freqency range of spectrum of all noise
##          frmin   = minimum freqency of good P, [0]: P, [1]: S
##          frmax   = maximum freqency of good P, [0]: P, [1]: S
##          goodP   = Bool variable for good P
##          goodS   = Bool variable for good S
    ## PARAMETERS FOR MULTI-TAPER
    nft=1024
    npi=3.0
    smlen = 11
    ## DETERMINE P SPECTRA
    for ii in range(pwindata.shape[0]):
        mtmresult=pymutt.mtft(pwindata[ii],dt=dt,npi=npi,
                              nwin=int(npi*2-1),paddedlen=nft)
        newspec=mtmresult['power'][1:]
        newfreq=(mtmresult['df']*np.arange(len(mtmresult['power'])))[1:]
        ## CONVERTING VELOCITY TO DISPLACEMENT BY DIVIDING BY 2*pi*f (Gubbins, p30)
        newspec=np.sqrt(newspec)/(2*np.pi*newfreq)
##        if smlen>0:
##            newspec=seis.smooth(newspec,smlen) ## SMOOTH THE SPECTRUM
        if ii==0:
            spec=newspec
            freq=newfreq
            finterv=mtmresult['df']
        else:
##            newspec=seis.smooth(newspec,21) ## SMOOTH THE NOISE SPECTRUM
            n_spec=newspec
            n_freq=newfreq
    ## DETERMINE S SPECTRA ON CHANNEL 1
    for ii in range(swindata1.shape[0]):
        mtmresult=pymutt.mtft(swindata1[ii],dt=dt,npi=npi,
                                nwin=int(npi*2-1),paddedlen=nft)
        newspec=mtmresult['power'][1:]
        newfreq=(mtmresult['df']*np.arange(len(mtmresult['power'])))[1:]
        newspec=np.sqrt(newspec)/(2*np.pi*newfreq)
##        if smlen>0:
##            newspec=seis.smooth(newspec,smlen) ## SMOOTH THE SPECTRUM
        if ii==0:   ## S WAVE
            spec=np.vstack((spec,newspec))
            freq=np.vstack((freq,newfreq))
            finterv=[finterv,mtmresult['df']]
        elif ii==1: ## S NOISE
##            newspec=seis.smooth(newspec,21) ## SMOOTH THE NOISE SPECTRUM
            n_spec=np.vstack((n_spec,newspec))
            n_freq=np.vstack((n_freq,newfreq))
        elif ii==2:  ## P CODA
            pcspec=newspec
    ## DETERMINE S SPECTRA ON CHANNEL 2
    for ii in range(swindata2.shape[0]):
        mtmresult=pymutt.mtft(swindata1[ii],dt=dt,npi=npi,
                                nwin=int(npi*2-1),paddedlen=nft)
        newspec=mtmresult['power'][1:]
        newfreq=(mtmresult['df']*np.arange(len(mtmresult['power'])))[1:]
        newspec=np.sqrt(newspec)/(2*np.pi*newfreq)
##        if smlen>0:
##            newspec=seis.smooth(newspec,smlen) ## SMOOTH THE SPECTRUM
        if ii==0:   ## S WAVE
            spec=np.vstack((spec,newspec))
            freq=np.vstack((freq,newfreq))
            finterv=[finterv,mtmresult['df']]
        elif ii==1: ## S NOISE
##            newspec=seis.smooth(newspec,21) ## SMOOTH THE NOISE SPECTRUM
            n_spec=np.vstack((n_spec,newspec))
            n_freq=np.vstack((n_freq,newfreq))
## USE P CODA
##        if ii==0:   ## S WAVE
##            sspec=newspec
##        elif ii==1: ## S NOISE
####            newspec=seis.smooth(newspec,21) ## SMOOTH THE NOISE SPECTRUM
##            n_spec=np.vstack((n_spec,newspec))
##            n_freq=np.vstack((n_freq,newfreq))
##        elif ii==2:  ## P CODA
####            pcspec=seis.smooth(newspec,21)  ## SMOOTH THE P CODA SPECTRUM
##            pcspec=newspec
##            spcspec=np.sqrt((sspec**2-pcspec**2).clip(min=0))
##            spec=np.vstack((spec,spcspec))
##            freq=np.vstack((freq,newfreq))
##            finterv=[finterv,mtmresult['df']]
    spec_px=spec
    freq_px=freq
    spec_sx=spec
    freq_sx=freq
    frmin=[6,6]
    frmax=[6,6]
    goodP=False
    goodS=False
    nsamp=[0,0]
    snrmed=[0,0]
    ## SINGAL-TO-NOISE RATIO
    snr=spec/n_spec
    if smlen>0:
        for ii in [0]:
            snr[ii]=seis.smooth(spec[ii],smlen)/seis.smooth(n_spec[ii],smlen)
##    lenspec=snr.shape[1]
    (begind,endind,frminp,frmaxp,frangep)=longseg(snr[0],snrcrtp,freq[0])
    frmin[0]=frminp
    frmax[0]=frmaxp
    if frangep<snrcrtp[1] or frminp>4:
        goodP=False
        goodS=False
        return spec_px,freq_px,spec_sx,freq_sx,spec,freq,n_spec,n_freq,frmin,frmax,goodP,goodS
    else:
        goodP=True
    spec_px=spec[0][begind:endind]
    freq_px=freq[0][begind:endind]
    spec_nx=n_spec[0][begind:endind]
##    spec_px=np.sqrt(spec_px**2-spec_nx**2)
    nsamp[0]=freq_px.shape[0]
    snr_px=snr[0][begind:endind]
    snrmed[0]=float(np.median(snr_px))
    coeffp=np.polyfit(freq_px,np.log(spec_px),1)
    synp=coeffp[1]+freq_px*coeffp[0]
##    residp=np.linalg.norm(np.log(synp)-np.log(spec_px))/np.sqrt(len(freq_px)-1)
    residp=seis.lincorrcoef(freq_px,np.log(spec_px))
    if coeffp[0]<0 and abs(residp)>=linresid[0]:
        goodP=True
    else:
        goodP=False
        goodS=False
        return spec_px,freq_px,spec_sx,freq_sx,spec,freq,n_spec,n_freq,frmin,frmax,goodP,goodS
    ## FIND THE LONGEST SEGMENT OF S SPECTRUM WITH SNR > SNRCRT
    if SDATA:
        (begind1,endind1,frmins1,frmaxs1,franges1)=longseg(snr[1],snrcrts,freq[1])
        (begind2,endind2,frmins2,frmaxs2,franges2)=longseg(snr[2],snrcrts,freq[2])
        if franges1>=franges2:
            begind =begind1
            endind =endind1
            frmins =frmins1
            frmaxs =frmaxs1
            franges=franges1
            sch=1
        else:
            begind =begind2
            endind =endind2
            frmins =frmins2
            frmaxs =frmaxs2
            franges=franges2
            sch=2
        frmin[1]=frmins
        frmax[1]=frmaxs
##        print(sta,frmins,frmaxs)
        if franges<snrcrts[1] or frmins>4:
            goodS=False
##            return spec_px,freq_px,spec_sx,freq_sx,spec,freq,n_spec,n_freq,frmin,frmax,goodP,goodS
        else:
            goodS=True
            spec_sx=spec[sch][begind:endind]
            freq_sx=freq[sch][begind:endind]
            spec_nx=n_spec[sch][begind:endind]
    ##        spec_sx=np.sqrt(spec_sx**2-spec_nx**2)
            nsamp[1]=freq_sx.shape[0]
            snr_sx=snr[sch][begind:endind]
            snrmed[1]=float(np.median(snr_sx))
            coeffs=np.polyfit(freq_sx,np.log(spec_sx),1)
            syns=coeffs[1]+freq_sx*coeffs[0]
    ##        resids=np.linalg.norm(np.log(syns)-np.log(spec_sx))/np.sqrt(len(freq_sx)-1)
            resids=seis.lincorrcoef(freq_sx,np.log(spec_sx))
            if coeffs[0]<0 and abs(resids)>=linresid[1]:
                goodS=True
            else:
                goodS=False
    if doplot:
        if SDATA and goodS:
            ## PLOT P AND S SIGNAL AND FITTING
            plt.figure(2)
            plt.clf()
            plt.subplot(2,2,1)
            plt.plot(freq[0],np.log(spec[0]),'b')   # P WAVE
            plt.plot(n_freq[0],np.log(n_spec[0]),'r')
            plt.plot(freq_px,np.log(spec_px),'k')
            plt.plot([frmin[0],frmin[0]],np.log([min(n_spec[0]),max(spec[0])]),'g')
            plt.plot([frmax[0],frmax[0]],np.log([min(n_spec[0]),max(spec[0])]),'g')
            plt.plot(freq_px,synp,'g--',linewidth=2)
##            plt.loglog(freq[0],spec[0],'b')   # P WAVE
##            plt.loglog(n_freq[0],n_spec[0],'r')
##            plt.loglog(freq_px,spec_px,'k')
##            plt.loglog([frmin[0],frmin[0]],[min(n_spec[0]),max(spec[0])],'g')
##            plt.loglog([frmax[0],frmax[0]],[min(n_spec[0]),max(spec[0])],'g')
##            plt.loglog(freq_px,synp,'g--',linewidth=2)
            plt.text(10,max(np.log(spec[0]))/2,'slope = %.4f' % coeffp[0])
            plt.text(10,max(np.log(spec[0]))/5,'residual = %.4f' % residp)
            plt.xlim([0.05,20])
            plt.ylabel('ln(Ap) on '+chan[2]+', nm/s')
            plt.title('Station: %s' % sta)            
            plt.subplot(2,2,2)
            plt.plot(freq[0],snr[0])
            plt.axhline(snrcrtp[0],color='g',linestyle='--')
            plt.axvline(frmin[0],color='g')
            plt.axvline(frmax[0],color='g')
            plt.xlim([0,20])
            plt.ylim([0,10])
            plt.ylabel('P Signal-to-Noise Ratio')
            plt.title('ORID = %d' % orid)            
            plt.subplot(2,2,3)
            plt.plot(freq[1],np.log(spec[1]),'b')   # S WAVE
##            plt.plot(freq[0],np.log(sspec),'y')   # S WAVE
##            plt.plot(freq[0],np.log(pcspec),'y--')   # P CODA
            plt.plot(n_freq[1],np.log(n_spec[1]),'r') # NOISE
########            plt.plot(freq2[1],np.log(spec2[1]),'g')
########            plt.plot(n_freq2[1],np.log(n_spec2[1]),'g--')
            plt.plot(freq_sx,np.log(spec_sx),'k')
            plt.plot([frmin[1],frmin[1]],np.log([min(n_spec[1]),max(spec[1])]),'g')
            plt.plot([frmax[1],frmax[1]],np.log([min(n_spec[1]),max(spec[1])]),'g')
            plt.plot(freq_sx,syns,'g--',linewidth=2)
##            plt.loglog(freq[1],spec[1],'b')   # S - P CODA
####            plt.loglog(freq[0],sspec,'y')   # S WAVE
##            plt.loglog(freq[0],pcspec,'y--')   # P CODA
##            plt.loglog(n_freq[1],n_spec[1],'r')
##            plt.loglog(freq_sx,spec_sx,'k')
##            plt.loglog([frmin[1],frmin[1]],[min(n_spec[1]),max(spec[1])],'g')
##            plt.loglog([frmax[1],frmax[1]],[min(n_spec[1]),max(spec[1])],'g')
##            plt.loglog(freq_sx,syns,'g--',linewidth=2)
            plt.xlim([0.05,20])
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('ln(As) on '+chan[sch]+', nm/s')
            plt.text(10,max(np.log(spec[1]))/2,'slope = %.4f' % coeffs[0])
            plt.text(10,max(np.log(spec[1]))/5,'residual = %.4f' % resids)
            plt.subplot(2,2,4)
            plt.plot(freq[1],snr[1])
            plt.axhline(snrcrts[0],color='g',linestyle='--')
            plt.axvline(frmin[1],color='g')
            plt.axvline(frmax[1],color='g')
            plt.ylim([0,10])
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('S Signal-to-Noise Ratio')
            plt.savefig(g.figdir+'/%d_%s_PSsnr.eps' % (orid,sta))
        else:
            ## PLOT ONLY P SIGNAL AND FITTING
            plt.figure(1)
            plt.clf()
            plt.subplot(2,2,1)
            plt.plot(freq[0],np.log(spec[0]),'b')
########        plt.plot(freq2[0],np.log(spec2[0]),'g')
########        plt.plot(n_freq2[0],np.log(n_spec2[0]),'g--')
            plt.plot(n_freq[0],np.log(n_spec[0]),'r')
            plt.plot(freq_px,np.log(spec_px),'k')
            plt.plot([frmin[0],frmin[0]],np.log([min(n_spec[0]),max(spec[0])]),'g')
            plt.plot([frmax[0],frmax[0]],np.log([min(n_spec[0]),max(spec[0])]),'g')
            plt.plot(freq_px,synp,'g--',linewidth=2)
            plt.text(10,max(np.log(spec[0]))/2,'slope = %.4f' % coeffp[0])
            plt.text(10,max(np.log(spec[0]))/10,'residual = %.4f' % residp)
########        plt.plot(matfreq_px,np.log(matspec_px),'k')
            plt.xlim([0,20])
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('ln(Ap) on '+chan[2]+', nm/s')
            plt.title('Station: %s' % sta)
            plt.subplot(2,2,2)
            plt.plot(freq[0],snr[0])
            plt.axhline(snrcrtp[0],color='g',linestyle='--')
            plt.axvline(frmin[0],color='g')
            plt.axvline(frmax[0],color='g')
            plt.xlim([0,20])
            plt.ylim([0,10])
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('P Signal-to-Noise Ratio')
            plt.title('ORID = %d' % orid)
            plt.subplot(2,2,3)
            plt.plot(freq[1],np.log(spec[1]),'b')   # S WAVE
            plt.plot(n_freq[1],np.log(n_spec[1]),'r') # NOISE
            plt.xlim([0.05,20])
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('ln(As) on '+chan[0]+', nm/s')
            plt.subplot(2,2,4)
            plt.plot(freq[1],snr[1])
            plt.axhline(snrcrts[0],color='g',linestyle='--')
            plt.axvline(frmin[1],color='g')
            plt.axvline(frmax[1],color='g')
            plt.ylim([0,10])
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('S Signal-to-Noise Ratio')
            plt.savefig(g.figdir+'/%d_%s_Psnr.eps' % (orid,sta))

    return spec_px,freq_px,spec_sx,freq_sx,spec,freq,n_spec,n_freq,frmin,frmax,goodP,goodS


def diffspec(pwindata,swindata1,swindata2,SDATA,orid,sta,snrcrt,linresid,chan,doplot,alpha):
## Calculate amplitude spectral ratio of S over P wave
## USAGE: (spec_px,freq_px,spec,freq,n_spec,n_freq,frmn,frmx,dtstar,
##         goodP2,goodS2)=tstarsub.diffspec(PWINDATA,SWINDATA1,SWINDATA2,dt,
##                   SDATA,orid,sta,snrcrtp,snrcrts,lincor,chan,doplot,alpha)
## INPUT:   pwindata  = P windowed data [0] and noise [1]
##          swindata1 = S windowed data [0] and noise [1] and P coda [2] on channel 1
##          swindata2 = S windowed data [0] and noise [1] and P coda [2] on channel 2
##          dt        = 1/sample rate
##          SDATA     = existence of S arrival (Bool variable)
##          orid      = origin id
##          sta       = station name
##          snrcrtp   = minimum [SNR,width,lowSNR] of good P in freqency domain
##          snrcrts   = minimum [SNR,width,lowSNR] of good S in freqency domain
##          lincor    = MINIMUM LINEAR CORRELATION COEFFICIENTS
##          chan      = CHANNELS
##          doplot    = 0: NO SNR PLOT, 1: ONLY S-P SNR PLOT, 2: P & S-P SNR PLOT
##          alpha     = Frequency dependence
## OUTPUT:  spec_px = spectrum of good P signal
##          freq_px = freqency range of spectrum of good P signal
##          spec    = spectrum of all signal
##          freq    = freqency range of spectrum of all signal
##          n_spec  = fspectrum of all noise
##          n_freq  = freqency range of spectrum of all noise
##          frmin   = minimum freqency of good P, [0]: P, [1]: S
##          frmax   = maximum freqency of good P, [0]: P, [1]: S
##          dtstar  = [dt*(S-P),std]
##          goodP   = Bool variable for good P
##          goodS   = Bool variable for good S
    ## PARAMETERS FOR MULTI-TAPER
    nft=1024
    npi=3.0
    dt = 0.025
    smlen = 11
    snrcrtp=[snrcrt[0],snrcrt[1],snrcrt[2]]
    snrcrts=[snrcrt[3],snrcrt[4],snrcrt[5]]
    ## DETERMINE P SPECTRA
    for ii in range(pwindata.shape[0]):
        mtmresult=pymutt.mtft(pwindata[ii],dt=dt,npi=npi,
                              nwin=int(npi*2-1),paddedlen=nft)
        newspec=mtmresult['power'][1:]
##        print(mtmresult['df'])
        newfreq=(mtmresult['df']*np.arange(len(mtmresult['power'])))[1:]
        ## CONVERTING VELOCITY TO DISPLACEMENT BY DIVIDING BY 2*pi*f (Gubbins, p30)
        newspec=np.sqrt(newspec)/(2*np.pi*newfreq)
##        if smlen>0:
##            newspec=seis.smooth(newspec,smlen) ## SMOOTH THE SPECTRUM
        if ii==0:
            spec=newspec
            freq=newfreq
            finterv=mtmresult['df']
        else:
            n_spec=newspec
            n_freq=newfreq
    ## DETERMINE S SPECTRA ON CHANNEL 1
    for ii in range(swindata1.shape[0]):
        mtmresult=pymutt.mtft(swindata1[ii],dt=dt,npi=npi,
                                nwin=int(npi*2-1),paddedlen=nft)
        newspec=mtmresult['power'][1:]
##        print(mtmresult['df'])
        newfreq=(mtmresult['df']*np.arange(len(mtmresult['power'])))[1:]
        newspec=np.sqrt(newspec)/(2*np.pi*newfreq)
##        if smlen>0:
##            newspec=seis.smooth(newspec,smlen) ## SMOOTH THE SPECTRUM
        if ii==0:   ## S WAVE
            spec=np.vstack((spec,newspec))
            freq=np.vstack((freq,newfreq))
            finterv=[finterv,mtmresult['df']]
        elif ii==1: ## S NOISE
            n_spec=np.vstack((n_spec,newspec))
            n_freq=np.vstack((n_freq,newfreq))
        elif ii==2:  ## P CODA
            pcspec=newspec
    ## DETERMINE S SPECTRA ON CHANNEL 2
    for ii in range(swindata2.shape[0]):
        mtmresult=pymutt.mtft(swindata1[ii],dt=dt,npi=npi,
                                nwin=int(npi*2-1),paddedlen=nft)
        newspec=mtmresult['power'][1:]
        newfreq=(mtmresult['df']*np.arange(len(mtmresult['power'])))[1:]
        newspec=np.sqrt(newspec)/(2*np.pi*newfreq)
##        if smlen>0:
##            newspec=seis.smooth(newspec,smlen) ## SMOOTH THE SPECTRUM
        if ii==0:   ## S WAVE
            spec=np.vstack((spec,newspec))
            freq=np.vstack((freq,newfreq))
            finterv=[finterv,mtmresult['df']]
        elif ii==1: ## S NOISE
            n_spec=np.vstack((n_spec,newspec))
            n_freq=np.vstack((n_freq,newfreq))
    spec_px=spec
    freq_px=freq
    spec_sx=spec
    freq_sx=freq
    spec_sp=spec
    freq_sp=freq
    frmin=[6,6]
    frmax=[6,6]
    goodP=False
    goodS=False
    nsamp=[0,0]
    snrmed=[0,0]
    dtstar=[0,0,0,0]
    ## SINGAL-TO-NOISE RATIO
    snr=spec/n_spec
    if smlen>0:
        for ii in [0]:
            snr[ii]=seis.smooth(spec[ii],smlen)/seis.smooth(n_spec[ii],smlen)
    (begind,endind,frminp,frmaxp,frangep)=longseg(snr[0],snrcrtp,freq[0])
    frmin[0]=frminp
    frmax[0]=frmaxp
    if frangep<snrcrtp[1] or frminp>4:
        goodP=False
        goodS=False
        return spec_px,freq_px,spec,freq,n_spec,n_freq,frmin,frmax,dtstar,goodP,goodS
    else:
        goodP=True
    spec_px=spec[0][begind:endind]
    freq_px=freq[0][begind:endind]
    spec_nx=n_spec[0][begind:endind]
##    spec_px=np.sqrt(spec_px**2-spec_nx**2)
    nsamp[0]=freq_px.shape[0]
    snr_px=snr[0][begind:endind]
    snrmed[0]=float(np.median(snr_px))
##    coeffp,residp = seis.polyodr(freq_px,np.log(spec_px),1)
##    residp=residp/np.sqrt(len(freq_px))
##    coeffp=np.polyfit(freq_px**(1-alpha),np.log(spec_px),1)
##    synp=coeffp[1]+(freq_px**(1-alpha))*coeffp[0]
##    tstarp=coeffp[0]/(-np.pi)
    coeffp=np.polyfit(freq_px,np.log(spec_px),1)
    synp=coeffp[1]+freq_px*coeffp[0]
##    residp=np.linalg.norm(np.log(synp)-np.log(spec_px))/np.sqrt(len(freq_px)-1)
    residp=seis.lincorrcoef(freq_px,np.log(spec_px))
    if coeffp[0]<0 and abs(residp)>=linresid[0]:
        goodP=True
    else:
        goodP=False
        goodS=False
        return spec_px,freq_px,spec,freq,n_spec,n_freq,frmin,frmax,dtstar,goodP,goodS
    ## FIND THE LONGEST SEGMENT OF S SPECTRUM WITH SNR > SNRCRT
    if SDATA:
        (begind1,endind1,frmins1,frmaxs1,franges1)=longseg(snr[1],snrcrts,freq[1])
        (begind2,endind2,frmins2,frmaxs2,franges2)=longseg(snr[2],snrcrts,freq[2])
        if franges1>=franges2:
            begind =begind1
            endind =endind1
            frmins =frmins1
            frmaxs =frmaxs1
            franges=franges1
            sch=0
        else:
            begind =begind2
            endind =endind2
            frmins =frmins2
            frmaxs =frmaxs2
            franges=franges2
            sch=1
        frmin[1]=frmins
        frmax[1]=frmaxs
        if franges<snrcrts[1] or frmins>4:
            goodS=False
##            return spec_px,freq_px,spec,freq,n_spec,n_freq,frmin,frmax,dtstar,goodP,goodS
        else:
            goodS=True
            spec_sx=spec[1][begind:endind]
            freq_sx=freq[1][begind:endind]
            spec_nx=n_spec[1][begind:endind]
    ##        spec_sx=np.sqrt(spec_sx**2-spec_nx**2)
            nsamp[1]=freq_sx.shape[0]
            snr_sx=snr[1][begind:endind]
            snrmed[1]=float(np.median(snr_sx))
    ##        coeffs,resids = seis.polyodr(freq_sx,np.log(spec_sx),1)
    ##        resids=resids/np.sqrt(len(freq_sx))
    ##        coeffs=np.polyfit(freq_sx**(1-alpha),np.log(spec_sx),1)
    ##        syns=coeffs[1]+(freq_sx**(1-alpha))*coeffs[0]
    ##        tstars=coeffs[0]/(-np.pi)
            coeffs=np.polyfit(freq_sx,np.log(spec_sx),1)
            syns=coeffs[1]+freq_sx*coeffs[0]
    ##        resids=np.linalg.norm(np.log(syns)-np.log(spec_sx))/np.sqrt(len(freq_sx)-1)
            resids=seis.lincorrcoef(freq_sx,np.log(spec_sx))
            if coeffs[0]>=0 or abs(resids)<linresid[1]:
                goodS=False
            else:
                goodS=True
##                return spec_px,freq_px,spec,freq,n_spec,n_freq,frmin,frmax,dtstar,goodP,goodS
                ## CALCULATE SPECTRAL RATIO OF S OVER P WAVES
                comindp=np.nonzero(np.in1d(freq_px,freq_sx))[0]
                cominds=np.nonzero(np.in1d(freq_sx,freq_px))[0]
                if len(comindp)>0:
                    spec_pcom=spec_px[comindp]
                    spec_scom=spec_sx[cominds]
                    freq_sp=freq_px[comindp]
                    spec_sp=spec_scom/spec_pcom
                    if (max(freq_sp)-min(freq_sp))>=linresid[3]:
##                        goodS=False
##                        return spec_px,freq_px,spec,freq,n_spec,n_freq,frmin,frmax,dtstar,goodP,goodS                
        ##            spec_sp=seis.smooth(spec_sp,21)
                        frmin[1]=min(freq_sp)
                        frmax[1]=max(freq_sp)
                        data=np.array([np.log(spec_sp)]).transpose()
                        G=np.array([-np.pi*freq_sp**(1-alpha)]).transpose()
                        G=np.hstack((np.ones((data.shape[0],1)),G))
                        model,residuS,mrank,sv=lstsq(G,data)
                        Ginv=np.linalg.inv(np.dot(G.transpose(),G))
                        synsp=np.dot(G,model)
                        dterr=np.sqrt(np.linalg.norm(data-synsp)**2/(data.shape[0]-2)
                                     *Ginv.diagonal()[1])  
                        dtstar=[float(model[1]),dterr,float(residuS**2/data.shape[0]),
                                float(residuS/data.shape[0])]
            ##            coeffsp=np.polyfit(freq_sp**(1-alpha),np.log(spec_sp),1)
            ####            coeffsp=np.polyfit(freq_sp,np.log(spec_sp),1)
            ##            synsp=coeffsp[1]+(freq_sp**(1-alpha))*coeffsp[0]
            ##            dtstar1=coeffsp[0]/(-np.pi)
            ##            residsp=np.linalg.norm(np.log(synsp)-np.log(spec_sp))/np.sqrt(len(freq_sp)-1)
                        residsp=seis.lincorrcoef(freq_sp,np.log(spec_sp))
##                        if dtstar[0]>0 and abs(residsp)>=linresid[2]:
                        if abs(residsp)>=linresid[2]:
                            goodS=True
                        else:
                            goodS=False
                    else:
                        goodS=False
        ##                    return spec_px,freq_px,spec,freq,n_spec,n_freq,frmin,frmax,dtstar,goodP,goodS
                else:
                    goodS=False
    ##                return spec_px,freq_px,spec,freq,n_spec,n_freq,frmin,frmax,dtstar,goodP,goodS
    if SDATA and goodS and doplot>0:
        ## PLOT P AND S SIGNAL AND FITTING
        plt.figure(2)
        plt.clf()
        plt.subplot(2,2,1)
        plt.plot(freq[0],np.log(spec[0]),'b')   # P WAVE
        plt.plot(n_freq[0],np.log(n_spec[0]),'r')
        plt.plot(freq_px,np.log(spec_px),'k')
        plt.plot([frmin[0],frmin[0]],np.log([min(n_spec[0]),max(spec[0])]),'g')
        plt.plot([frmax[0],frmax[0]],np.log([min(n_spec[0]),max(spec[0])]),'g')
        plt.plot(freq_px,synp,'g--',linewidth=2)
        plt.text(10,max(np.log(spec[0]))/2,'Slope = %.4f' % coeffp[0])
        plt.text(10,max(np.log(spec[0]))/5,'lincor = %.4f' % residp)
        plt.xlim([0.05,20])
        plt.ylabel('ln(Ap) on '+chan[2]+', nm/s')
        plt.title('Station: %s' % sta)            
        plt.subplot(2,2,2)
        plt.plot(freq[0],snr[0])
        plt.axhline(snrcrtp[0],color='g',linestyle='--')
        plt.axvline(frmin[0],color='g')
        plt.axvline(frmax[0],color='g')
        plt.xlim([0,20])
        plt.ylim([0,10])
        plt.ylabel('P Signal-to-Noise Ratio')
        plt.title('ORID = %d' % orid)            
        plt.subplot(2,2,3)
        plt.plot(freq[1],np.log(spec[1]),'b')   # S WAVE
        plt.plot(n_freq[1],np.log(n_spec[1]),'r') # NOISE
        plt.plot(freq_sx,np.log(spec_sx),'k')
        plt.plot([frmins,frmins],np.log([min(n_spec[1]),max(spec[1])]),'g')
        plt.plot([frmaxs,frmaxs],np.log([min(n_spec[1]),max(spec[1])]),'g')
        plt.plot(freq_sx,syns,'g--',linewidth=2)
        plt.xlim([0.05,20])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('ln(As) on '+chan[sch]+', nm/s')
        plt.text(10,max(np.log(spec[1]))/2,'Slope = %.4f' % coeffs[0])
        plt.text(10,max(np.log(spec[1]))/5,'lincor = %.4f' % resids)
        plt.subplot(2,2,4)
        plt.plot(freq_sp,np.log(spec_sp),'b')   # S-P RATIO
        plt.plot([frmin[1],frmin[1]],np.log([min(spec_sp),max(spec_sp)]),'g')
        plt.plot([frmax[1],frmax[1]],np.log([min(spec_sp),max(spec_sp)]),'g')
        plt.plot(freq_sp,synsp,'g--',linewidth=2)
        plt.xlim([0.05,((np.floor(frmax[1]/5)+1)*5)])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('ln(As)/ln(Ap)')
        plt.text(frmin[1]+1,0.75*max(np.log(spec_sp))+
                 0.25*min(np.log(spec_sp)),'dt* = %.4f' % dtstar[0])
        plt.text(frmin[1]+1,(max(np.log(spec_sp))+
                min(np.log(spec_sp)))/2,'std = %.4f' % dterr)
##            plt.text(frmin[1]+1,0.75*max(np.log(spec_sp))+
##                     0.25*min(np.log(spec_sp)),'dt* = %.4f' % dtstar1)
        plt.text(frmin[1]+1,0.25*max(np.log(spec_sp))+
                0.75*min(np.log(spec_sp)),'lincor = %.4f' % residsp)
        plt.savefig(g.figdir+'/%d_%s_PSsnr.eps' % (orid,sta))
    elif goodP and doplot==2:
        ## PLOT ONLY P SIGNAL AND FITTING
        plt.figure(1)
        plt.clf()
        plt.subplot(2,2,1)
        plt.plot(freq[0],np.log(spec[0]),'b')
########        plt.plot(freq2[0],np.log(spec2[0]),'g')
########        plt.plot(n_freq2[0],np.log(n_spec2[0]),'g--')
        plt.plot(n_freq[0],np.log(n_spec[0]),'r')
        plt.plot(freq_px,np.log(spec_px),'k')
        plt.plot([frmin[0],frmin[0]],np.log([min(n_spec[0]),max(spec[0])]),'g')
        plt.plot([frmax[0],frmax[0]],np.log([min(n_spec[0]),max(spec[0])]),'g')
        plt.plot(freq_px,synp,'g--',linewidth=2)
        plt.text(10,max(np.log(spec[0]))/2,'Slope = %.4f' % coeffp[0])
        plt.text(10,max(np.log(spec[0]))/10,'lincor = %.4f' % residp)
########        plt.plot(matfreq_px,np.log(matspec_px),'k')
        plt.xlim([0,20])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('ln(Ap) on '+chan[2]+', nm/s')
        plt.title('Station: %s' % sta)
        plt.subplot(2,2,2)
        plt.plot(freq[0],snr[0])
        plt.axhline(snrcrtp[0],color='g',linestyle='--')
        plt.axvline(frmin[0],color='g')
        plt.axvline(frmax[0],color='g')
        plt.xlim([0,20])
        plt.ylim([0,10])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('P Signal-to-Noise Ratio')
        plt.title('ORID = %d' % orid)
        plt.subplot(2,2,3)
        plt.plot(freq[1],np.log(spec[1]),'b')   # S WAVE
        plt.plot(n_freq[1],np.log(n_spec[1]),'r') # NOISE
        plt.xlim([0.05,20])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('ln(As) on '+chan[0]+', nm/s')
        plt.subplot(2,2,4)
        plt.plot(freq[1],snr[1])
        plt.ylim([0,10])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('S Signal-to-Noise Ratio')
        plt.savefig(g.figdir+'/%d_%s_Psnr.eps' % (orid,sta))
    ## IN FUTURE MAY NEED TO ADD P/S PART
    return spec_px,freq_px,spec,freq,n_spec,n_freq,frmin,frmax,dtstar,goodP,goodS


def plotspec(saving,sta,orid,POS,lnmomen,fc,alpha,icase,sitedata=0):
## PLOT AMPLITUDE SPECTRUM
    if POS.upper()=='P':
        ind=0
        xmax=10
        ymin=-5
        ymax=10
        textx=6
    elif POS.upper()=='S':
        ind=1
        xmax=4
        ymin=4
        ymax=12
        textx=2.5
    else:
        raise ValueError, "P or S wave?"
    corr=saving['corr'][ind]
    spec=saving['spec'][ind]
    freq=saving['freq'][ind]
    n_spec=saving['nspec'][ind]
    n_freq=saving['nfreq'][ind]
    frmin=saving[2]['frmin'][ind]
    frmax=saving[2]['frmax'][ind]
    invtstar=saving[icase]['tstar'][ind]
    synspec=(corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar)/(1+(freq/fc)**2))
    if POS.upper()=='S':
        invtstarP=saving[icase]['tstar'][0]
        ttP=saving['Ptt']
        ttS=saving['Stt']
        QpQs=2.25
        invtstar2=invtstarP*QpQs*ttS/ttP
        synspec2=(corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar2)/(1+(freq/fc)**2))
        QpQs=1.75
        invtstar2=invtstarP*QpQs*ttS/ttP
        synspec3=(corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar2)/(1+(freq/fc)**2))
    indx=np.all([(freq>=frmin),(freq<frmax)],axis=0)
    specx=spec[indx]
    freqx=freq[indx]
    synx=synspec[indx]
    resid=(1-((np.linalg.norm(np.log(synx)-np.log(specx)))**2/(len(freqx)-1)
             /np.var(np.log(specx))))    
    df=abs(freq[1]-freq[0])
    nlowf=0
    narea=0
    for ifreq in range(len(freq)):
        if (freq[ifreq]>frmax and freq[ifreq]<15):
            if (np.log(synspec[ifreq])<np.log(spec[ifreq]) or
                 np.log(synspec[ifreq])>np.log(spec[ifreq])+1):
                narea=narea+np.log(spec[ifreq])-np.log(synspec[ifreq])
            if np.log(synspec[ifreq])>np.log(spec[ifreq])+2:
                nlowf=nlowf+5
            elif np.log(synspec[ifreq])>np.log(spec[ifreq])+1:
                nlowf=nlowf+1
    if narea<-10 and nlowf*df>3:
        resid=0
##    if POS.upper()=='P':
##        est=saving[icase]['est'][ind].transpose()[0]-np.log(1+(freqx/fc)**2)+np.log(corr)
##        dat=saving[icase]['dat'][ind].transpose()[0]-np.log(1+(freqx/fc)**2)+np.log(corr)
##    elif POS.upper()=='S':
##        est=saving[icase]['est'][ind].transpose()[0]-np.log(1+(freqx/fc)**2)+np.log(corr)+lnmomen
##        dat=saving[icase]['dat'][ind].transpose()[0]-np.log(1+(freqx/fc)**2)+np.log(corr)+lnmomen
    plt.figure(99)
    plt.clf()
    plt.plot(freq,np.log(spec),'b')
    if not isinstance(sitedata, int):
        plt.plot(freq,np.log(spec)-sitedata,'b--')
    plt.plot(n_freq,np.log(n_spec),'r')
    plt.plot([frmin,frmin],np.log([min(n_spec),max(spec)]),'g')
    plt.plot([frmax,frmax],np.log([min(n_spec),max(spec)]),'g')
    plt.plot(freq,np.log(synspec),'g--',linewidth=2)
    if POS.upper()=='S':
        plt.plot(freq,np.log(synspec2),'k:',linewidth=2)
        plt.plot(freq,np.log(synspec3),'k--',linewidth=2)
##    plt.loglog(freq,spec,'b')
##    plt.loglog(n_freq,n_spec,'r')
##    plt.loglog([frmin,frmin],[min(n_spec),max(spec)],'g')
##    plt.loglog([frmax,frmax],[min(n_spec),max(spec)],'g')
##    plt.loglog(freq,synspec,'g--',linewidth=2)
##    plt.plot(freqx,est,'k--')
##    plt.plot(freqx,dat,'k')
    plt.text(textx,max(np.log(spec))/2,'t* = %.2f' % invtstar)
##    plt.text(10,max(np.log(spec))/4,'curve fitting = %.4f' % resid)
##    plt.text(10,max(np.log(spec))/7,'misfit = %.4f' % saving[icase]['misfit'][ind])
##    plt.text(10,max(np.log(spec))/2.5,'sd = %.4f' % saving[icase]['err'][ind])
    plt.xlim([0,xmax])
#    plt.ylim(np.log([min(n_spec),max(spec)]))
    plt.ylim([ymin,ymax])
#    plt.ylim([-5,np.log(max(spec))])
##    plt.ylim([10**(int(np.log(min(min(spec),min(n_spec))))-2),
##                   10**(int(np.log(max(max(spec),max(n_spec))))+1)])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('ln(A%s), nm/s' % POS.lower())
    plt.title('ORID = %d  Station: %s' % (orid,sta))
    plt.savefig(g.figdir+'/%d_%s_%sspectrum.pdf' % (orid,sta,POS.upper()))
    return

def fitting(saving,sta,orid,POS,lnmomen,fc,alpha,icase):
## CALCULATE HOW WELL THE SYNTHETIC SPECTRUM FITS THE DATA
## IF THE FITTING CURVE IS BELOW THE NOISE, THEN resid = 999999.
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError, "P or S wave?"
    corr=saving['corr'][ind]
    spec=saving['spec'][ind]
    freq=saving['freq'][ind]
    frmin=saving[2]['frmin'][ind]
    frmax=saving[2]['frmax'][ind]
    invtstar=saving[icase]['tstar'][ind]
    synspec=(corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar)/(1+(freq/fc)**2))
    indx=np.all([(freq>=frmin),(freq<frmax)],axis=0)
    specx=spec[indx]
    freqx=freq[indx]
    synx=synspec[indx]
##    resid=np.linalg.norm(np.log(synx)-np.log(specx))/(len(freqx))
    resid=(1-((np.linalg.norm(np.log(synx)-np.log(specx)))**2/(len(freqx)-1)
             /np.var(np.log(specx))))    
    df=abs(freq[1]-freq[0])
    nlowf=0
    narea=0
    for ifreq in range(len(freq)):
        if (freq[ifreq]>frmax and freq[ifreq]<15):
            if (np.log(synspec[ifreq])<np.log(spec[ifreq]) or
                 np.log(synspec[ifreq])>np.log(spec[ifreq])+1):
                narea=narea+np.log(spec[ifreq])-np.log(synspec[ifreq])
            if np.log(synspec[ifreq])>np.log(spec[ifreq])+2:
                nlowf=nlowf+5
            elif np.log(synspec[ifreq])>np.log(spec[ifreq])+1:
                nlowf=nlowf+1
    if narea<-10 and nlowf*df>3:
        resid=0
    return resid

def calresspec(saving,sta,orid,POS,lnmomen,fc,alpha):
## CALCULATE RESIDUAL SPECTRUM FOR EACH STATION
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError, "P or S wave?"
    freq_x = saving[2][POS.lower()][0]
    spec_x = saving[2][POS.lower()][1]
    correc = saving['corr'][ind]
    invtstar=saving[3]['tstar'][ind]
    righthand = lnmomen-np.pi*freq_x*(freq_x**(-alpha)*invtstar)
    resspec = np.array([np.log(spec_x)-np.log(correc)
                  +np.log(1+(freq_x/fc)**2)-righthand])
    resratio = resspec/righthand*100
    resspec = np.vstack((freq_x,resspec))
    resspec = np.vstack((resspec,resratio))
    resspec = resspec.transpose()
#     print(resspec)
#     exit()
    return resspec

def buildd(saving,stalst,fc,POS,icase,lnM=0):
## Build data matrix
##      d = [ln(A1)-ln(C1)+ln(1+(f1i/fc)**2),                   ##
##           ln(A2)-ln(C2)+ln(1+(f2i/fc)**2),                   ##
##           ln(AM)-ln(CM)+ln(1+(fMi/fc)**2)]                   ##
## INPUT:   saving - saved spectrum for each station: saving[sta][1]['p']
##          stalst - list of used stations
##          fc     - corner frequency
##          POS    - 'P' or 'S'
##          icase  - 1: high quality for finding best fc and alpha
##                   2: low quality for t* inversion
##                   3: low quality for t* inversion without bad fitting
##          lnM    - when POS='S', log of seismic moment
## OUTPUT:  data - data matrix for t* inversion
    if icase==3:
        icase=2
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError, "P or S wave?"
    for ista in range(len(stalst)):
        sta=stalst[ista]
        freq_x = saving[sta][icase][POS.lower()][0]
        spec_x = saving[sta][icase][POS.lower()][1]
        correc = saving[sta]['corr'][ind]
        stad = np.array([np.log(spec_x)-np.log(correc)
                        +np.log(1+(freq_x/fc)**2)-lnM]).transpose()
##        print(sta,POS,max(stad),lnM)
        if ista==0:
            data=stad
        else:
            data=np.vstack((data,stad))
    return data

def buildd_site(saving,stalst,fc,allsite,POS,icase,lnM=0):
## Build data matrix with site correction
##      d = [ln(A1)-ln(C1)+ln(1+(f1i/fc)**2),                   ##
##           ln(A2)-ln(C2)+ln(1+(f2i/fc)**2),                   ##
##           ln(AM)-ln(CM)+ln(1+(fMi/fc)**2)]                   ##
## INPUT:   saving - saved spectrum for each station: saving[sta][1]['p']
##          stalst - list of used stations
##          fc     - corner frequency
##          allsite - site effects of all stations
##          POS    - 'P' or 'S'
##          icase  - 1: high quality for finding best fc and alpha
##                   2: low quality for t* inversion
##                   3: low quality for t* inversion without bad fitting
##          lnM    - when POS='S', log of seismic moment
## OUTPUT:  data - data matrix for t* inversion
    if icase==3:
        icase=2
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError, "P or S wave?"
    for ista in range(len(stalst)):
        sta=stalst[ista]
        freq_x = saving[sta][icase][POS.lower()][0]
        spec_x = saving[sta][icase][POS.lower()][1]
        correc = saving[sta]['corr'][ind]
        if sta in allsite.keys():
            sitedata=allsite[sta]
            sitefreq=sitedata[:,0]
            sitespec=sitedata[:,1]
#             print(sitefreq,freq_x)
#             ind1=np.nonzero(freq_x[0]==sitefreq)[0][0]
            ind1=np.nonzero(abs(freq_x[0]-sitefreq)<0.01)[0][0]
            ind2=ind1+len(freq_x)
            sitespec_x=sitespec[ind1:ind2]
            stad = np.array([np.log(spec_x)-np.log(correc)
                        +np.log(1+(freq_x/fc)**2)-sitespec_x-lnM]).transpose()
        else:
            stad = np.array([np.log(spec_x)-np.log(correc)
                        +np.log(1+(freq_x/fc)**2)-lnM]).transpose()
##        print(sta,POS,max(stad),lnM)
        if ista==0:
            data=stad
        else:
            data=np.vstack((data,stad))
    return data
    
def buildG(saving,stalst,alpha,POS,icase):
## Build G matrix
##      G = [[1, -pi*f1i*f1i**(-alpha), 0, ..., 0],             ##
##           [1, 0, -pi*f2i*f2i**(-alpha), ..., 0],             ##
##           [1, 0, 0, ..., -pi*fMi*fMi**(-alpha)]]             ##
## INPUT:   saving - saved frequency for each station: saving[sta][1]['p']
##          stalst - list of used stations
##          alpha  - alpha value(s)
##          POS    - 'P' or 'S'
##          icase  - 1: high quality for finding best fc and alpha
##                   2: low quality for t* inversion
##                   3: low quality for t* inversion without bad fitting
## OUTPUT:  G - G matrix for t* inversion
    if icase==3:
        icase=2
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError, "P or S wave?"
    for ista in range(len(stalst)):
        sta=stalst[ista]
        for alco in range(len(alpha)):
            freq_x = saving[sta][icase][POS.lower()][0]
            exponent = -1*np.pi*freq_x*(freq_x**(-alpha[alco]))
            exponent = np.array([exponent]).transpose()
            if alco==0:
                Gblock = np.atleast_3d(exponent)
            else:
                Gblock = np.dstack((Gpblock,exponent))
        if ista==0:
            G=Gblock
        else:
            oldblock=np.hstack((G,np.zeros((G.shape[0],1,len(alpha)))))
            newblock=np.hstack((np.zeros((Gblock.shape[0],G.shape[1],len(alpha))),Gblock))
            G=np.vstack((oldblock,newblock))
    if POS.upper()=='P':
        G=np.hstack((np.ones((G.shape[0],1,len(alpha))),G))
    return G

def buildGmos(saving,stalst,alpha,POS,icase):
## Build G matrix
##      G = [[1, -pi*f1i*f1i**(-alpha), 0, ..., 0],             ##
##           [1, 0, -pi*f2i*f2i**(-alpha), ..., 0],             ##
##           [1, 0, 0, ..., -pi*fMi*fMi**(-alpha)]]             ##
## INPUT:   saving - saved frequency for each station: saving[sta][1]['p']
##          stalst - list of used stations
##          alpha  - alpha value(s)
##          POS    - 'P' or 'S'
##          icase  - 1: high quality for finding best fc and alpha
##                   2: low quality for t* inversion
##                   3: low quality for t* inversion without bad fitting
## OUTPUT:  G - G matrix for t* inversion
    if icase==3:
        icase=2
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError, "P or S wave?"
    for ista in range(len(stalst)):
        sta=stalst[ista]
        for alco in range(len(alpha)):
            freq_x = saving[sta][icase][POS.lower()][0]
            exponent = -1*np.pi*freq_x*(freq_x**(-alpha[alco]))
            exponent = np.array([exponent]).transpose()
            if alco==0:
                Gblock = np.atleast_3d(exponent)
            else:
                Gblock = np.dstack((Gpblock,exponent))
        if ista==0:
            G=Gblock
        else:
            oldblock=np.hstack((G,np.zeros((G.shape[0],1,len(alpha)))))
            newblock=np.hstack((np.zeros((Gblock.shape[0],G.shape[1],len(alpha))),Gblock))
            G=np.vstack((oldblock,newblock))
    G=np.hstack((np.ones((G.shape[0],1,len(alpha))),G))
    return G
