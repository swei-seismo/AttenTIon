#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##Sub functions for t* inversion 
import numpy as np
import os
from scipy.signal import *
from scipy.linalg import lstsq
import matplotlib.pyplot as plt
import matplotlib as mpl
import tstar_parameters as tp
import multitaper as mtm
import seissign as seis

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.size']=8
mpl.rcParams['axes.formatter.limits']=[-2,2]

def readseismo(pretime,dt,orid,sta,chan):
    """
    Read seismograms in txt format from './data/processedSeismograms/'.
    USAGE: (dd,time,flag) = readseismo(pretime,dt,orid,sta,chan)
    INPUT:  (pretime)  ==> Seconds before arrival to subset
            (dt)       ==> Scalar: 1/(samplerate of data)
            (subdb)    ==> Subset of database
            (orid)     ==> Origin ID of event
            (sta)      ==> Station name
            (chan)     ==> Channel names (i.e. ['BHE','BHN','BHZ'])
    OUTPUT: (dd)       ==> (3Xn) numpy.ndarray: [(E_CHAN),(N_CHAN),(Z_CHAN)]
            (time)     ==> (1Xn) numpy.ndarray: (Relative Time)
            (flag)     ==> If reading successful
    """
    for ichan in range(len(chan)):
        sacfl=tp.sacdir+'/%s/%s_%s.%s.txt'%(orid, orid.split('_')[1],sta, chan[ichan])
        if not os.path.isfile(sacfl):
            print('ERROR: %s does not exist' % sacfl)
            dd=np.array(range(18000))
            tt=np.array(range(18000))
            flag=False
            return dd,tt,flag
        ddchan=np.fromstring(("".join(open(sacfl).readlines()[30:])),sep=' ')
        nn=ddchan.size
        #yurong 06/26/22: if error occurs here, elements in nn (1x3) are possibly inconsistent,
        #need to find the smallest nn instead of the original size
        if ichan==0:
            dd=ddchan
        else:
            dd=np.vstack((dd,ddchan))
    flag=True
    tt=-pretime+dt*np.array(range(nn)) 
    
    return dd,tt,flag

def fixwin(dd,tt,param,chan,ARRIV,orid,sta):
    """
    Window seismic data in time domain.
    USAGE: (p_dd,s_dd,pn_dd,sn_dd,NOS2)
            =tstarsub.fixwin(dd,tt,dt,ARRIV[ii+1],prewin,WLP,WLS,NOS,doplot,orid,sta)
    INPUT: (dd)      ==> (3Xn) numpy.ndarray: [(E_CHAN),(N_CHAN),(Z_CHAN)]
        (tt)      ==> (1Xn) numpy.ndarray: (Relative Time)
        (dt)      ==> Scalar: 1/(samplerate of data)
        (chan)    ==> List of channal names
        (ARRIV)   ==> Arrival Dictionary with Keys (see getdbinfo)
        (prewin)  ==> Seconds before P[0]/S[1] arrival for windowing
        (WLP)     ==> Window lenght for P in seconds
        (WLS)     ==> Window lenght for S in seconds
        (SDATA)   ==> Existence of S arrival
        (doplot)  ==> Bool variable for plotting spectrum
        (orid)    ==> origin id
        (sta)     ==> station name
    OUTPUT: (p_dd)   ==> Windowed P data, starting prewin seconds before P arrival
            (s_dd)   ==> Windowed S data
            (pn_dd)  ==> P noise data, WL length ending prewin seconds before P arrival
            (sn_dd)  ==> S noise data
            (pd_dd)  ==> P coda data, right before S arrivel

    """
    ## WINDOWING P
    pch=chan.index(ARRIV['pchan'])
    pind=np.all([(tt>=-1*param['prewin'][0]),(tt<=(param['WLP']-param['prewin'][0]))],axis=0)
    p_dd=detrend(dd[pch][pind]-np.mean(dd[pch][pind]))
    pnind=np.all([(tt<=-1*param['prewin'][0]),(tt>=(-param['WLP']-param['prewin'][0]))],axis=0)
    pn_dd=detrend(dd[pch][pnind]-np.mean(dd[pch][pnind]))
    ## MAKING SURE P WAVE AND NOISE HAVE SAME SIZE
    if p_dd.size>pn_dd.size:
        p_dd=p_dd[0:(pn_dd.size)]
    elif p_dd.size<pn_dd.size:
        pn_dd=p_dd[0:(p_dd.size)]   

    ## Plot P window data
    if param['doplotseis']:
        scalefac=1e3
        plt.figure(3)
        plt.clf()
        plt.subplot(2,1,2)
        plt.xlabel(ARRIV['pchan']+' (seconds)')
        plt.ylabel('Velocity Amplitude (nm/s)')
        plt.title('Station: %s' % sta)
        # plt.plot(tt,dd[pch])
        tmpdd=dd[pch]
        b, a = butter(4, 0.01, btype='highpass')
        filtmpdd = filtfilt(b, a, tmpdd)
        plt.plot(tt,filtmpdd)
        # plt.xlim([-param['WLP']-param['prewin'][0]-3,param['WLP']-param['prewin'][0]+3])
        plt.xlim([-10,10])
        plt.ylim([-7e4,7e4])
        # plt.ylim([np.floor(min(filtmpdd)/scalefac)*scalefac,np.ceil(max(filtmpdd)/scalefac)*scalefac])
        plt.axvline(-1*param['prewin'][0],color='g')
        plt.axvline(param['WLP']-param['prewin'][0],color='g')
        plt.axvline(-param['WLP']-param['prewin'][0],color='g')
        plt.text(-param['WLP']-param['prewin'][0],min(p_dd),'Noise', fontsize = 15)
        plt.text(-1*param['prewin'][0],min(p_dd),'P wave',fontsize = 15)
        plt.savefig(tp.figdir1+'/%s_%s_Pwindata.pdf' % (orid,sta))  

    ## WINDOWING S ON BOTH HORIZONTAL CHANNELS
    # if SDATA:
    #     sminusp=ARRIV['stime']-ARRIV['ptime']
    #     if sminusp<(param['WLS']+param['WLP']+param['prewin'][0]-param['prewin'][1]+param['gaps']):
    #         print('P & S arrivels are too close - proceeding as if no S pick')
    #         SDATA=False 
        #yurong: possible error source for S wave: sminusp
        # if round(sminusp/0.025,5) == int(sminusp/0.025):
        #     sminusp += 0.001
    ## yurong: work for S data
    return p_dd,pn_dd

def longseg(snr,snrcrtpara,freq,minf=0.05,maxf=15.0):
    """
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
    """
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
    # print(frmin,frmax,frange,m)
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
                # print(frmin2,frmax2,frange2,snrcrtpara[1])
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

def dospec(pwindata, orid, sta, param, chan, data_quality):
    """
    Calculate amplitude spectrum of windowed waveform using multi-taper method
    USAGE: (spec_px,freq_px,spec_sx,freq_sx,spec,freq,n_spec,n_freq,frmn,frmx,
            goodP1,goodS1)=tstarsub.dospec(PWINDATA,SWINDATA1,SWINDATA2,dt,
                        SDATA,orid,sta,snrcrtp,snrcrts,lincor,chan,doplot)
    INPUT:  pwindata  = P windowed data [0] and noise [1]
            swindata1 = S windowed data [0] and noise [1] and P coda [2] on channel 1
            swindata2 = S windowed data [0] and noise [1] and P coda [2] on channel 2
            dt        = 1/sample rate
            SDATA     = existence of S arrival (Bool variable)
            orid      = origin id
            sta       = station name
            snrcrtp   = minimum [SNR,width,lowSNR] of good P in freqency domain
            snrcrts   = minimum [SNR,width,lowSNR] of good S in freqency domain
            lincor    = MINIMUM LINEAR CORRELATION COEFFICIENTS
            chan      = CHANNELS
            doplot    = Bool variable for plotting spectrum
    OUTPUT:  spec_px = spectrum of good P signal
            freq_px = freqency range of spectrum of good P signal
            spec_sx = spectrum of good S signal
            freq_sx = freqency range of spectrum of good S signal
            spec    = spectrum of all signal
            freq    = freqency range of spectrum of all signal
            n_spec  = fspectrum of all noise
            n_freq  = freqency range of spectrum of all noise
            frmin   = minimum freqency of good P, [0]: P, [1]: S
            frmax   = maximum freqency of good P, [0]: P, [1]: S
            goodP   = Bool variable for good P
            goodS   = Bool variable for good S
    """
    ## PARAMETERS FOR MULTI-TAPER
    nft=1024
    npi=3.0
    smlen = 11
    spec_px, freq_px = [], []
    ## DETERMINE P SPECTRA

    for ii in range(pwindata.shape[0]):
        mtmresult=mtm.mtspec(pwindata[ii],param['dt'],time_bandwidth=5, nfft = nft)
        newspec=mtmresult[0][1:]
        newfreq=mtmresult[1][1:]
        ## CONVERTING VELOCITY TO DISPLACEMENT BY DIVIDING BY 2*pi*f (Gubbins, p30)
        newspec=np.sqrt(newspec)/(2*np.pi*newfreq)
        if ii==0:
            pspec=newspec
            pfreq=newfreq
        else:
            pn_spec=newspec
            pn_freq=newfreq
    ##yurong: S wave multitaper

    frmin, frmax = [6,6], [6,6]
    goodP, goodS = False, False
    psnr=seis.smooth(pspec,smlen)/seis.smooth(pn_spec,smlen)
    (begind,endind,frminp,frmaxp,frangep)=longseg(psnr,param['snrcrtp'+str(data_quality)],pfreq)
    frmin[0], frmax[0] = frminp, frmaxp
    if frangep<param['snrcrtp'+str(data_quality)][1] or frminp>param['frminp']:
        goodP=False
        goodS=False
        return goodP, goodS, spec_px, freq_px, pspec, pn_spec, pfreq, pn_freq, frmin, frmax
    else:
        goodP=True
    spec_px=pspec[begind:endind]
    freq_px=pfreq[begind:endind]
    coeffp=np.polyfit(freq_px,np.log(spec_px),1)
    synp=coeffp[1]+freq_px*coeffp[0]
    residp=seis.lincorrcoef(freq_px,np.log(spec_px))

    if coeffp[0]<0 and abs(residp)>=param['lincor'][0]:
        goodP=True
    else:
        goodP=False
        goodS=False
        return goodP, goodS, spec_px, freq_px, pspec, pn_spec, pfreq, pn_freq, frmin, frmax
    if data_quality != 1 and param['doplotsnr'] == True: ##plot seismic signals in frequency domain
        ##yurong: plot S wave signal
        ## PLOT ONLY P SIGNAL AND FITTING
            plt.figure(1)
            plt.clf()
            plt.subplot(2,2,1)
            plt.plot(pfreq,np.log(pspec),'b')
            # plt.plot(freq2[0],np.log(spec2[0]),'g')
            # plt.plot(n_freq2[0],np.log(n_spec2[0]),'g--')
            plt.plot(pn_freq,np.log(pn_spec),'r')
            plt.plot(freq_px,np.log(spec_px),'k')
            plt.plot([frmin[0],frmin[0]],np.log([min(pn_spec),max(pspec)]),'g')
            plt.plot([frmax[0],frmax[0]],np.log([min(pn_spec),max(pspec)]),'g')
            plt.plot(freq_px,synp,'g--',linewidth=2)
            plt.text(10,max(np.log(pspec))/2,'slope = %.4f' % coeffp[0])
            plt.text(10,max(np.log(pspec))/10,'lincorr = %.4f' % residp)
            plt.xlim([0,20])
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('ln(Ap) on '+chan[2]+', nm/s')
            plt.title('Station: %s' % sta)
            plt.subplot(2,2,2)
            plt.plot(pfreq,psnr)
            plt.axhline(param['snrcrtp'+str(data_quality)][0],color='g',linestyle='--')
            plt.axvline(frmin[0],color='g')
            plt.axvline(frmax[0],color='g')
            plt.xlim([0,20])
            plt.ylim([0,250])
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('P Signal-to-Noise Ratio')
            plt.title('%s' % orid)
            # plt.subplot(2,2,3)
            # if sch == 1:
            #     plt.plot(s1freq,np.log(s1spec),'b')   # S WAVE
            #     plt.plot(s1n_freq,np.log(s1n_spec),'r') # NOISE
            # elif sch == 2:
            #     plt.plot(s2freq,np.log(s2spec),'b')   # S WAVE
            #     plt.plot(s2n_freq,np.log(s2n_spec),'r') # NOISE
            # plt.xlim([0.05,20])
            # plt.xlabel('Frequency (Hz)')
            # plt.ylabel('ln(As) on '+chan[sch-1]+', nm/s')
            # plt.subplot(2,2,4)
            # if sch == 1:
            #     plt.plot(s1freq,s1snr)
            # elif sch == 2:
            #     plt.plot(s2freq,s2snr)
            # plt.axhline(snrcrts[0],color='g',linestyle='--')
            # plt.axvline(frmin[1],color='g')
            # plt.axvline(frmax[1],color='g')
            # plt.ylim([0,30])
            # plt.xlabel('Frequency (Hz)')
            # plt.ylabel('S Signal-to-Noise Ratio')
            plt.savefig(tp.figdir2+'/%s_%s_Psnr.pdf' % (orid,sta))



    return goodP, goodS, spec_px, freq_px, pspec, pn_spec, pfreq, pn_freq, frmin, frmax

def buildd(saving,stalst,ORIG,POS,icase,source_para,fc,allsite,add_site,lnM=0):
    """
    Build data matrix
        d = [ln(A1)-ln(C1)+ln(1+(f1i/fc)**2),                   ##
            ln(A2)-ln(C2)+ln(1+(f2i/fc)**2),                   ##
            ln(AM)-ln(CM)+ln(1+(fMi/fc)**2)]                   ##
    INPUT:   saving - saved spectrum for each station: saving[sta][1]['p']
            stalst - list of used stations
            fc     - corner frequency
            allsite - site effects of all stations
            POS    - 'P' or 'S'
            icase  - 
                    1: high quality for finding best fc and alpha
                    2: low quality for t* inversion
                    3: low quality for t* inversion without bad fitting
            lnM    - when POS='S', log of seismic moment
    OUTPUT:  data - data matrix for t* inversion
    """
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError("P or S wave?")
    
    mw, mo = ORIG['mw'], ORIG['mo']
    if source_para == 2:
        lnM = np.log(mo)

    for ista in range(len(stalst)):
        sta=stalst[ista]
        freq_x = saving[sta][icase][POS.lower()][0]
        spec_x = saving[sta][icase][POS.lower()][1]
        correc = saving[sta]['corr'][ind]
        if add_site == 0:
            sitespec_x = 0
        elif sta in allsite.keys():
            sitedata=allsite[sta]
            sitefreq=sitedata[:,0]
            sitespec=sitedata[:,1]
            ind1=np.nonzero(abs(freq_x[0]-sitefreq)<0.01)[0][0]
            ind2=ind1+len(freq_x)
            sitespec_x=sitespec[ind1:ind2]
        else: ## 1)add_site = 1 and 2)sta is not in allsite.leys()
            sitespec_x = 0
        stad = np.array([np.log(spec_x)-np.log(correc)
                        +np.log(1+(freq_x/fc)**2)-sitespec_x-lnM]).transpose()
        # print(sta,POS,max(stad),lnM)
        if ista==0:
            data=stad
        else:
            data=np.vstack((data,stad))
    return data

def buildG(saving,stalst,alpha,POS,icase,source_para):
    """
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
    """    
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError("P or S wave?")
    for ista in range(len(stalst)):
        sta=stalst[ista]
        freq_x = saving[sta][icase][POS.lower()][0]
        exponent = -1*np.pi*freq_x*(freq_x**(-alpha))
        exponent = np.array([exponent]).transpose()
        Gblock = np.atleast_3d(exponent)

        if ista==0:
            G=Gblock
        else:
            oldblock=np.hstack((G,np.zeros((G.shape[0],1,1))))
            newblock=np.hstack((np.zeros((Gblock.shape[0],G.shape[1],1)),Gblock))
            G=np.vstack((oldblock,newblock))
    if source_para == 1:    ## grid search for moment magnitude
        if POS.upper()=='P':
            G=np.hstack((np.ones((G.shape[0],1,1)),G))

    return G

def fitting(saving,sta,ORIG,POS,alpha,lnmomen,icase):
## CALCULATE HOW WELL THE SYNTHETIC SPECTRUM FITS THE DATA
## IF THE FITTING CURVE IS BELOW THE NOISE, THEN resid = 999999.
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError("P or S wave?")
    corr=saving[sta]['corr'][ind]
    freq, spec=saving[sta][icase][POS.lower()][0], saving[sta][icase][POS.lower()][1]
    invtstar=saving[sta][icase]['tstar'][ind]
    synx=(corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar)/(1+(freq/ORIG['fc'])**2))

    resid=(1-((np.linalg.norm(np.log(synx)-np.log(spec)))**2/(len(freq)-1)/np.var(np.log(spec)))) 
    ##    RESID = 1-√(∑((ln(A(synthetic))-ln(A(observed)))^2))/(n-1)/σ(ln(A(observed))) Bill Menke, 'geophysical data analysis discrete inverse theory'
    resid1=(np.linalg.norm(np.log(synx)-np.log(spec))/(len(freq)))/0.30
    ##    L2-NORM MISFIT = √(∑((ln(A(synthetic))-ln(A(observed)))^2))/n. Here 0.30 is just for normalization
    resid2=(np.linalg.norm(np.log(synx)-np.log(spec),ord=1)/(len(freq)))/1.50
    ##    L1-NORM MISFIT = (∑|ln(A(synthetic))-ln(A(observed))|)/n. Here 1.50 is just for normalization
    resid3=(1-2*np.sum(np.log(synx)*np.log(spec))/(np.sum((np.log(spec))**2)+np.sum((np.log(synx))**2)))/0.8
    ##    CORRELATIVE FUNCTION. 1-2*∑(ln(A(synthetic))*ln(A(observed)))/∑((ln(A(synthetic))^2)+(ln(A(observed)))^2). Here 0.80 is just for normalization

    ##yurong: what do these lines mean?
    # for ifreq in range(len(freq)):
    #     if (freq[ifreq]>frmax and freq[ifreq]<15):
    #         if (np.log(synspec[ifreq])<np.log(spec[ifreq]) or
    #              np.log(synspec[ifreq])>np.log(spec[ifreq])+1):
    #             narea=narea+np.log(spec[ifreq])-np.log(synspec[ifreq])
    #         if np.log(synspec[ifreq])>np.log(spec[ifreq])+2:
    #             nlowf=nlowf+5
    #         elif np.log(synspec[ifreq])>np.log(spec[ifreq])+1:
    #             nlowf=nlowf+1
    # if narea<-10 and nlowf*df>3:
    #     resid=0
    return resid

def calresspec(saving,POS,lnmomen,fc,alpha):
## CALCULATE RESIDUAL SPECTRUM FOR EACH STATION
    if POS.upper()=='P':
        ind=0
    elif POS.upper()=='S':
        ind=1
    else:
        raise ValueError("P or S wave?")
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
    return resspec

def plotspec(saving,sta,orid,POS,lnmomen,fc,alpha,icase,sitedata=0):
    """PLOT AMPLITUDE SPECTRUM
    """
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
        raise ValueError("P or S wave?")
    corr=saving['corr'][ind]
    spec=saving['spec'][ind]
    freq=saving['freq'][ind]
    n_spec=saving['nspec'][ind]
    n_freq=saving['nfreq'][ind]
    frmin=saving[2]['frmin'][ind]
    frmax=saving[2]['frmax'][ind]
    invtstar=saving[icase]['tstar'][ind]
    synspec=(corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar)/(1+(freq/fc)**2))
    ##yurong: S wave work
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
    ## yurong: the same question in fitting function
    # nlowf=0
    # narea=0
    # for ifreq in range(len(freq)):
    #     if (freq[ifreq]>frmax and freq[ifreq]<15):
    #         if (np.log(synspec[ifreq])<np.log(spec[ifreq]) or
    #              np.log(synspec[ifreq])>np.log(spec[ifreq])+1):
    #             narea=narea+np.log(spec[ifreq])-np.log(synspec[ifreq])
    #         if np.log(synspec[ifreq])>np.log(spec[ifreq])+2:
    #             nlowf=nlowf+5
    #         elif np.log(synspec[ifreq])>np.log(spec[ifreq])+1:
    #             nlowf=nlowf+1
    # if narea<-10 and nlowf*df>3:
    #     resid=0
    plt.figure(1)
    plt.clf()
    if not isinstance(sitedata, int):
        plt.plot(freq,np.log(spec)-sitedata,'b--')
    plt.plot(freq,np.log(spec),'b')
    plt.plot(n_freq,np.log(n_spec),'r')
    plt.plot(freq,np.log(synspec),'g--',linewidth=2)
    plt.plot([frmin,frmin],np.log([min(n_spec),max(spec)]),'g')
    plt.plot([frmax,frmax],np.log([min(n_spec),max(spec)]),'g')
    
    ## yurong: S wave work
    # if POS.upper()=='S':
    #     plt.plot(freq,np.log(synspec2),'k:',linewidth=2)
    #     plt.plot(freq,np.log(synspec3),'k--',linewidth=2)
    plt.text(textx,max(np.log(spec))/2,'t* = %.2f' % invtstar)
    plt.xlim([0,xmax])
    plt.ylim([ymin,ymax])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('ln(A%s), nm/s' % POS.lower())
    plt.title('ORID = %s  Station: %s' % (orid,sta))
    plt.savefig(tp.figdir3+'/%s_%s_%sspectrum.pdf' % (orid,sta,POS.upper()))
    return