import os
import glob
import numpy as np
from scipy.linalg import lstsq
from scipy.signal import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import multitaper as mtm
import seissign as seis

def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def get_event_info(orid, catalog):
##  INPUT:  (orid)        ==> Event ID
##          (catalog)     ==> Path to the earthquake catalog
##  OUTPUT: (event_dic)   ==> dictionary contains event info with keys:
##              [orid] = origin id (orid)
##              [time] = origin time (unix seconds)
##              [lat]  = hypocenter latitude
##              [lon]  = hypocenter longitude
##              [dep]  = hypocenter depth in km
    event_dic = {}
    command = "cat %s | awk '/%s/ {print}'" %(catalog, orid)
    line = os.popen(command).read().strip()
    event_dic["orid"]  = orid
    event_dic["time"] = line.split()[0]
    event_dic["lat"] = float(line.split()[2])
    event_dic["lon"] = float(line.split()[1])
    event_dic["dep"] = float(line.split()[3])
    
    return event_dic


def get_station_info(orid,database_path):
##  INPUT:  (orid)          ==> Event ID
##          (database_path) ==> Path to the database
##  OUTPUT: (station_dic)   ==> dictionary contains station info:
    station_dic = {}
    database = database_path+"/"+"database_%s" %(orid)
    fdatabase = open(database,"r")
    for line in fdatabase.readlines()[1:]:
        line = line.strip()
        line_sp = line.split()
        sta = line_sp[1]
        station_dic[sta] = {}
        # saclst KNETWK KSTNM T0 T1 DIST BAZ GCARC DELTA B f *.SAC (create database)
        station_dic[sta]["net"] = line_sp[0]          # network name
        station_dic[sta]["sta"] = line_sp[1]          # station name
        station_dic[sta]["ptime"] = float(line_sp[2]) # P arrival time
        station_dic[sta]["stime"] = float(line_sp[3]) # S arrival time
        station_dic[sta]["dist"] = float(line_sp[4])  # distance in km
        station_dic[sta]["baz"] = float(line_sp[5])   # backazimuth
        station_dic[sta]["gcarc"] = float(line_sp[6]) # distance in degree
        station_dic[sta]["dt"] = float(line_sp[7])    # sampling rate
        station_dic[sta]["b"] = float(line_sp[8])     # begin time 
    fdatabase.close()
    
    return station_dic


def get_channels(sacdir, orid, net, sta):
##  INPUT:  (sacdir)  ==>  Path to Processed seismograms
##          (orid)    ==>  Event ID
##          (net)     ==>  Network name
##          (sta)     ==>  Station name
##  OUTPUT: (chan)    ==>  A list contains channel names
    sacfiles = sacdir+"/%s/%s.%s.*.SAC" %(orid, net, sta)
    command = "saclst KCMPNM f %s | awk -F ' ' '{print $2}'" %(sacfiles)
    output = os.popen(command).read().strip().split('\n')
    chan = sorted(output) # confirm chan[0]-E, chan[1]-N,chan[2]-Z
    
    return chan 


def readseismo(sacdir, pretime, dt, orid, net, sta, chan):
##  Read seismograms
##  INPUT:  (sacdir)   ==>  Path to Processed seismograms
##          (pretime)  ==>  Seconds before the P wave arrival
##          (dt)       ==>  Sampling rate
##          (orid)     ==>  Event ID
##          (net)      ==>  Network name
##          (sta)      ==>  Station name
##          (chan)     ==>  Channel names (like ["BHE","BHN","BHZ"])
##  OUTPUT: (dd)       ==>  (3Xn) numpy.ndarray: dd[0]-E, dd[1]-N, dd[2]-Z
##          (tt)       ==>  (1Xn) numpy.ndarray: (relative time)
##          (flag)     ==>  (If reading successful)
    for i in range(len(chan)):
        sacfile = glob.glob("%s/%s/%s.%s*%s.txt" %(sacdir,orid,net,sta,chan[i]))[0]
        if not os.path.isfile(sacfile):
            print('ERROR: %s does not exist' % sacfile)
            dd = np.array(range(18000))
            tt = np.array(range(18000))
            flag = False
            return dd, tt, flag
        else:
            ddchan = np.fromstring(("".join(open(sacfile).readlines()[30:])),sep=' ')
            nn = ddchan.size
            if i == 0:
                nmin = nn
            else:
                if nn < nmin:
                    nmin = nn
    
    for ichan in range(len(chan)):
        sacfile = glob.glob("%s/%s/%s.%s*%s.txt" %(sacdir,orid,net,sta,chan[i]))[0]
        ddchan = np.fromstring(("".join(open(sacfile).readlines()[30:])),sep=' ')
        ddchan = ddchan[:nmin]
        if ichan == 0:
            dd = ddchan
        else:
            dd = np.vstack((dd,ddchan))
    flag = True
    tt = -pretime+dt*np.array(range(nmin))

    return dd, tt, flag


def fixwin(dd,tt,dt,chan,ARRIV,prewin,WLP,WLS,SDATA,doplotseis,orid,net,sta,figdir1,PStooclose):
##  Windows seismic data
##  USAGE:  (p_dd,pn_dd,s_dd1,sn_dd1,pc_dd1,s_dd2,sn_dd2,pc_dd2,SDATA,PStooclose) = tstar.sub
##          (dd,tt,dt,chan,ARRIV,prewin,WLP,WLS,SDATA,doplotseis,orid,net,sta,figdir1,PStooclose)
##  INPUT:  (dd)         ==>  (3Xn) numpy.ndarray: dd[0]-E, dd[1]-N, dd[2]-Z
##          (tt)         ==>  (1Xn) numpy.ndarray: (relative time)
##          (dt)         ==>  Sampling rate
##          (chan)       ==>  A list contains channel names (E/N/Z)
##          (ARRIV)      ==>  Arrival dictionary with keys
##          (prewin)     ==>  Seconds before P[0]/S[1] arrival for windowing
##          (WLP)        ==>  Window length for P in seconds
##          (WLS)        ==>  Window length for S in seconds
##          (SDATA)      ==>  Existence of S arrival (bool variable) 
##          (doplotseis) ==>  Bool variable for plotting seismogram
##          (orid)       ==>  Event ID
##          (net)        ==>  Network name
##          (sta)        ==>  Station Name
##          (figdir1)    ==>  Path to plotting seismograms
##          (PStooclose) ==>  Number of P and S arrivals are too close
##  OUTPUT: (p_dd)       ==>  P signal data
##          (pn_dd)      ==>  P noise data
##          (s_dd1)      ==>  S signal data on channel 1
##          (sn_dd1)     ==>  S noise data on channel 1
##          (pc_dd1)     ==>  P coda data on channel 1
##          (s_dd2)      ==>  S signal data on channel 2
##          (sn_dd2)     ==>  S noise data on channel 2
##          (pc_dd2)     ==>  P coda data on channel 2
##          (SDATA)      ==>  Existence of S arrival (bool variable) 
##          (PStooclose) ==>  Number of P and S arrivals are too close
    gaps = 0.5 # seconds between S noise and signal

    ## Window P wave
    pch = chan.index(ARRIV['pchan'])
    pind = np.all([(tt >= -1*prewin[0]),(tt <= (WLP-prewin[0]))],axis=0)
    p_dd = detrend(dd[pch][pind]-np.mean(dd[pch][pind]))
    pnind = np.all([(tt <= -1*prewin[0]),(tt >= (-WLP-prewin[0]))],axis=0)
    pn_dd = detrend(dd[pch][pnind]-np.mean(dd[pch][pnind]))
    ## Make sure P signal and noise have same size
    if p_dd.size > pn_dd.size:
        p_dd = p_dd[0:(pn_dd.size)]
    elif p_dd.size < pn_dd.size:
        pn_dd = pn_dd[0:(p_dd.size)]
    ## Plot P wave seismogram
    if doplotseis:
        seis_figure_dir = figdir1+"/%s" %(orid)
        mkdir(seis_figure_dir)
        plt.figure(1)
        plt.clf()
        plt.figure(figsize=(12,6),dpi=200)
        plt.plot(tt,dd[pch],"blue",linewidth=1.5)
        plt.xlim([-10,10])
        plt.ylim(-1.2*max(abs(p_dd)),1.2*max(abs(p_dd)))
        plt.xlabel("Time relative to P arrival time (s)")
        plt.ylabel("Velocity Amplitude (nm/s)")
        plt.axvline(-1*prewin[0],color="green",linewidth=1.0)
        plt.axvline(WLP-prewin[0],color="green",linewidth=1.0)
        plt.axvline(-WLP-prewin[0],color="green",linewidth=1.0)
        plt.text(-WLP-prewin[0],-1.1*max(abs(p_dd)),"Noise", fontsize=15)
        plt.text(-1*prewin[0],-1.1*max(abs(p_dd)),"P wave",fontsize=15)
        figure_title = "%s %s.%s.%s" %(orid,net,sta,chan[pch])
        plt.title(figure_title)
        plt.savefig(seis_figure_dir+'/%s_%s_Pwindata.pdf' % (net,sta))
    
    ## Window S wave on both horizontal channels
    if SDATA:
        sminusp = ARRIV['stime']-ARRIV['ptime']
        if sminusp < (WLS+WLP+prewin[0]-prewin[1]+gaps):
            PStooclose += 1
            print('P & S arrivals are too close - proceeding as if no S pick')
            SDATA = False
        if round(sminusp/dt,5) == int(sminusp/dt):
            sminusp += 0.001
    if SDATA:
        ## Horizontal channel 1
        sch = 0
        sind = np.all([(tt>=(sminusp-prewin[1])),(tt<=(sminusp+WLS-prewin[1]))],axis=0)
        s_dd1 = detrend(dd[sch][sind]-np.mean(dd[sch][sind]))
        ## Noise is defined as gaps s before S arrival
        snind = np.all([(tt<=(sminusp-prewin[1]-gaps)),(tt>=(sminusp-WLS-prewin[1]-gaps))],axis=0)
        sn_dd1 = detrend(dd[sch][snind]-np.mean(dd[sch][snind]))
        ## P code is defined as right before S arrival
        pcind = np.all([(tt<=(sminusp-prewin[1])),(tt>=(sminusp-WLS-prewin[1]))],axis=0)
        pc_dd1 = detrend(dd[sch][pcind]-np.mean(dd[sch][pcind]))
        ## Make sure S wave, noise and P code wave have same size
        minlen = min(s_dd1.size, sn_dd1.size, pc_dd1.size)
        s_dd1 = s_dd1[0:minlen]
        sn_dd1 = sn_dd1[0:minlen]
        pc_dd1 = pc_dd1[0:minlen]
        ## Horizontal channel 2
        sch = 1
        sind = np.all([(tt>=(sminusp-prewin[1])),(tt<=(sminusp+WLS-prewin[1]))],axis=0)
        s_dd2 = detrend(dd[sch][sind]-np.mean(dd[sch][sind]))
        ## Noise is defined as gaps s before S arrival
        snind = np.all([(tt<=(sminusp-prewin[1]-gaps)),(tt>=(sminusp-WLS-prewin[1]-gaps))],axis=0)
        sn_dd2 = detrend(dd[sch][snind]-np.mean(dd[sch][snind]))
        ## P code is defined as right before S arrival
        pcind = np.all([(tt<=(sminusp-prewin[1])),(tt>=(sminusp-WLS-prewin[1]))],axis=0)
        pc_dd2 = detrend(dd[sch][pcind]-np.mean(dd[sch][pcind]))
        ## Make sure S wave, noise and P code wave have same size
        minlen = min(s_dd2.size,sn_dd2.size,pc_dd2.size)
        s_dd2 = s_dd2[0:minlen]
        sn_dd2 = sn_dd2[0:minlen]
        pc_dd2 = pc_dd2[0:minlen]
        ## Plot S wave seismogram
        if doplotseis:
            ## Plot seismogram for horizontal channel 1
            plt.figure(2)
            plt.clf()
            plt.figure(figsize=(12,6),dpi=200)
            plt.subplot(2,1,1)
            plt.plot(tt,dd[0],"blue",linewidth=1.5)
            plt.xlim([sminusp-10,sminusp+10])
            plt.ylim(-1.2*max(abs(s_dd1)),1.2*max(abs(s_dd1)))
            plt.axvline(sminusp-prewin[1],color="green")
            plt.axvline(sminusp+WLS-prewin[1],color="green")
            plt.axvline(sminusp-prewin[1]-gaps,color="green")
            plt.axvline(sminusp-WLS-prewin[1]-gaps,color="green")
            plt.axvline(-1*prewin[1],color="green")
            plt.text(-prewin[1],min(s_dd1),"P arrival")
            plt.text(sminusp-WLS-prewin[1]-gaps,-1.1*max(abs(s_dd1)),"Noise",fontsize=15)
            plt.text(sminusp-prewin[1],-1.1*max(abs(s_dd1)),"S wave",fontsize=15)
            figure_title = "%s %s.%s.%s" %(orid,net,sta,chan[0])
            plt.title(figure_title)
            ## Plot seismogram for horizontal channel 2
            plt.subplot(2,1,2)
            plt.plot(tt,dd[1],"blue",linewidth=1.5)
            plt.xlim([sminusp-10,sminusp+10])
            plt.ylim(-1.2*max(abs(s_dd2)),1.2*max(abs(s_dd2)))
            plt.axvline(sminusp-prewin[1],color="green")
            plt.axvline(sminusp+WLS-prewin[1],color="green")
            plt.axvline(sminusp-prewin[1]-gaps,color="green")
            plt.axvline(sminusp-WLS-prewin[1]-gaps,color="green")
            plt.axvline(-1*prewin[1],color="green")
            plt.text(-prewin[1],min(s_dd2),"P arrival")
            plt.text(sminusp-WLS-prewin[1]-gaps,-1.1*max(abs(s_dd2)),"Noise",fontsize=15)
            plt.text(sminusp-prewin[1],-1.1*max(abs(s_dd2)),"S wave",fontsize=15)
            figure_title = "%s %s.%s.%s" %(orid,net,sta,chan[1])
            plt.title(figure_title)
            plt.savefig(seis_figure_dir+'/%s_%s_Swindata.pdf' % (net,sta))
    else:
        s_dd1 = p_dd
        sn_dd1 = pn_dd
        pc_dd1 = p_dd
        s_dd2 = p_dd
        sn_dd2 = pn_dd
        pc_dd2 = p_dd
    
    return p_dd, pn_dd, s_dd1, sn_dd1, pc_dd1, s_dd2, sn_dd2, pc_dd2, SDATA, PStooclose


def longseg(snr, snrcrtpara, freq): 
##  FIND THE LONGEST SEGMENT OF SPECTRA WITH SNR > SNRCRT
##  USAGE: (begind,endind,frmin,frmax,frange) = tstarsub.longseg(snr,snrcrtp,freq)
##  INPUT:  (snr)        ==>  SIGNAL-TO-NOISE RATIO
##          (snrcrtpara) ==>  [MINIMUM SIGNAL-TO-NOISE, MINIMUM LENGTH OF SEGMENT]
##          (freq)       ==>  FREQUENCY
##  OUTPUT: (begind)     ==>  INDEX OF THE BEGIN OF THE LONGEST SEGMENT
##          (endind)     ==>  INDEX OF THE END OF THE LONGEST SEGMENT
##          (frmin)      ==>  MINIMUM FREQUENCY OF THE LONGEST SEGMENT
##          (frmax)      ==>  MAXIMUM FREQUENCY OF THE LONGEST SEGMENT
##          (frange)     ==>  frmax - frmin
    ## TAKE SPECTRUM < maxf (50 Hz) and > minf (0.05 Hz)
    maxf = 50
    minf = 0.35
    lenspec = len([ifreq for ifreq in freq if (ifreq<maxf and ifreq>=minf)])
    ind1 = int(min(np.nonzero(freq>=minf)[0]))
    ind2 = int(max(np.nonzero(freq<maxf)[0]))
    w = 0
    m = []
    bindex = []
    eindex = []
    snrcrt = snrcrtpara[0] # set SNR
    for kk in range(ind1+1,lenspec):
        if snr[kk]<snrcrt and snr[kk-1]>=snrcrt and kk==1:        # only first > crt
            w = 1
            m.append(w)
            bindex.append(kk-w)
            eindex.append(kk-1)
            w = 0
        elif snr[kk]>=snrcrt and snr[kk-1]>=snrcrt and kk==1:     # at first and continuously > crt
            w = w+2
        elif snr[kk]>=snrcrt and snr[kk-1]<snrcrt and kk>=1 and kk<(lenspec-1):    # begin of continuously > crt
            w = w+1
        elif snr[kk]>=snrcrt and snr[kk-1]>=snrcrt and kk>1 and kk<(lenspec-1):   # continuously >= crt
            w = w+1
        elif snr[kk]<snrcrt and snr[kk-1]>=snrcrt and kk>1 and kk<=(lenspec-1):    # end of continuously > crt
            m.append(w)
            bindex.append(kk-w)
            eindex.append(kk-1)
            w = 0
        elif snr[kk]<snrcrt and snr[kk-1]<snrcrt and kk>=1 and kk<=(lenspec-1):     # continuously < crt
            w = 0
        elif snr[kk]>=snrcrt and snr[kk]>=snrcrt and kk==(lenspec-1):     # at last and continuously > crt
            w = w+1
            m.append(w)
            bindex.append(kk-w+1)
            eindex.append(kk)
        elif snr[kk]>=snrcrt and snr[kk]<snrcrt and kk==(lenspec-1):      # only last > crt
            w = 1
            m.append(w)
            bindex.append(kk-w+1)
            eindex.append(kk)
    if len(m) == 0:
        frange = 0
        frmin = 6
        frmax = 6
        begind = 0
        endind = 0
        return begind,endind,frmin,frmax,frange
    
    ## FIND THE LONGEST SEGMENT
    longest = m.index(max(m))
    frmin = freq[bindex[longest]]
    frmax = freq[eindex[longest]]
    frange = frmax-frmin

    begind = bindex[longest]
    endind = eindex[longest]

    ## EXTEND FREQUENCY BAND TO lowSNR
    if snrcrtpara[2]<snrcrtpara[0]:
        if begind>ind1+1:
             while snr[begind-1]<snr[begind] and snr[begind-1]>snrcrtpara[2] and begind-1>ind1+1:
                 begind=begind-1
        if endind<ind2-1:
            while snr[endind+1]<snr[endind] and snr[endind+1]>snrcrtpara[2] and endind+1<ind2-1:
                endind = endind+1
        frmin = freq[begind]
        frmax = freq[endind]
        frange = frmax-frmin

    
    return begind,endind,frmin,frmax,frange

def dospec(pwindata,swindata1,swindata2,dt,SDATA,orid,net,sta,snrcrtp,snrcrts,
           linresid,chan,doplotsnr,figdir2):
## Calculate amplitude spectrum of windowed waveform using multi-taper method
## USAGE: (spec_px,freq_px,spec_sx,freq_sx,spec,freq,n_spec,n_freq,frmin,frmax,
##         goodP1,goodS1) = tstarsub.dospec(PWINDATA,SWINDATA1,SWINDATA2,dt,
##                         SDATA,orid,net,sta,snrcrtp,snrcrts,lincor,chan,doplotsnr)
## INPUT:  (pwindata)  ==> P windowed data [0] and noise [1]
##         (swindata1) ==> S windowed data [0] and noise [1] and P coda [2] on channel 1
##         (swindata2) ==> S windowed data [0] and noise [1] and P coda [2] on channel 2
##         (dt)        ==>  1/sample rate
##         (SDATA)     ==>  existence of S arrival (Bool variable)
##         (orid)      ==>  origin id
##         (net)       ==>  network name
##         (sta)       ==>  station name
##         (snrcrtp)   ==>  minimum [SNR,width,lowSNR] of good P in freqency domain
##         (snrcrts)   ==>  minimum [SNR,width,lowSNR] of good S in freqency domain
##         (lincor)    ==>  MINIMUM LINEAR CORRELATION COEFFICIENTS
##         (chan)      ==>  CHANNELS
##         (doplotsnr) ==>  Bool variable for plotting spectrum
##         (figdir2)   ==>  Path to plotting SNR of spectra
## OUTPUT: (spec_px)   ==>  spectrum of good P signal
##         (freq_px)   ==>  freqency range of spectrum of good P signal
##         (spec_sx)   ==>  spectrum of good S signal
##         (freq_sx)   ==>  freqency range of spectrum of good S signal
##         (spec)      ==>  spectrum of all signal
##         (freq)      ==>  freqency range of spectrum of all signal
##         (n_spec)    ==>  fspectrum of all noise
##         (n_freq)    ==>  freqency range of spectrum of all noise
##         (frmin)     ==>  minimum freqency of good P, [0]: P, [1]: S
##         (frmax)     ==>  maximum freqency of good P, [0]: P, [1]: S
##         (goodP)     ==>  Bool variable for good P
##         (goodS)     ==>  Bool variable for good S

    smlen = 11
    residp = 100
    resids = 100
    pspec = []
    s1spec = []
    s2spec = []
    pn_spec = []
    s1n_spec = []
    s2n_spec = []
    pfreq = []
    s1freq = []
    s2freq = []
    pn_freq = []
    s1n_freq = []
    s2n_freq = []
    sch = -1

    ## DETERMINE P SPECTRA
    for ii in range(pwindata.shape[0]):
# =============================================================================
        mtmresult = mtm.sine_psd(pwindata[ii],dt)
# =============================================================================
        newspec = mtmresult[0][1:]
        newfreq = mtmresult[1][1:]
        ## Converting velocity to displacement by dividing by 2*pi*f (Gubbins, p30)
        newspec = np.sqrt(newspec)/(2*np.pi*newfreq)
        if ii == 0:
            pspec = newspec
            pfreq = newfreq
            finterv = mtmresult[1]
        else:
            pn_spec = newspec
            pn_freq = newfreq

    ## DETERMINE S SPECTRA ON CHANNEL 1
    for ii in range(swindata1.shape[0]):
# =============================================================================
        mtmresult = mtm.sine_psd(swindata1[ii],dt)
# =============================================================================
        newspec = mtmresult[0][1:]
        newfreq = mtmresult[1][1:]
        newspec = np.sqrt(newspec)/(2*np.pi*newfreq)
        if ii == 0:   ## S WAVE
            s1spec = newspec
            s1freq = newfreq
        elif ii == 1: ## S NOISE
            s1n_spec = newspec
            s1n_freq = newfreq
        elif ii == 2:  ## P CODA
            pcspec = newspec
    
    ## DETERMINE S SPECTRA ON CHANNEL 2
    for ii in range(swindata2.shape[0]):
# =============================================================================
        mtmresult = mtm.sine_psd(swindata2[ii],dt)
# =============================================================================
        newspec = mtmresult[0][1:]
        newfreq = mtmresult[1][1:]
        newspec = np.sqrt(newspec)/(2*np.pi*newfreq)
        if ii == 0:   ## S WAVE
            s2spec = newspec
            s2freq = newfreq
        elif ii == 1: ## S NOISE
            s2n_spec = newspec
            s2n_freq = newfreq

    spec_px = pspec
    freq_px = pfreq
    spec_sx = s1spec
    freq_sx = s1freq
    frmin = [6,6]
    frmax = [6,6]
    goodP = False
    goodS = False
    nsamp = [0,0]
    snrmed = [0,0] # Median snr
    
    ## SINGAL-TO-NOISE RATIO
    psnr = pspec/pn_spec
    s1snr = s1spec/s1n_spec
    s2snr = s2spec/s2n_spec
    ## FIND THE LONGEST SEGMENT OF P SPECTRUM WITH SNR > SNRCRT
    if smlen > 0:
        psnr = seis.smooth(pspec,smlen)/seis.smooth(pn_spec,smlen)
    (begind,endind,frminp,frmaxp,frangep) = longseg(psnr,snrcrtp,pfreq)
    frmin[0] = frminp
    frmax[0] = frmaxp
    if frangep < snrcrtp[1] or frminp > 4:
        goodP = False
        goodS = False
        return sch,residp,resids,spec_px,freq_px,spec_sx,freq_sx,pspec, \
               s1spec,s2spec,pfreq,s1freq,s2freq,pn_spec,s1n_spec, \
               s2n_spec,pn_freq,s1n_freq,s2n_freq,frmin,frmax,goodP,goodS
    else:
        goodP = True
    spec_px = pspec[begind:endind]
    freq_px = pfreq[begind:endind]
    spec_nx = pn_spec[begind:endind]
    nsamp[0] = freq_px.shape[0]
    snr_px = psnr[begind:endind]
    snrmed[0] = float(np.median(snr_px))
    coeffp = np.polyfit(freq_px,np.log(spec_px),1) # linear fitting
    synp = coeffp[1]+freq_px*coeffp[0]
    residp = seis.lincorrcoef(freq_px,np.log(spec_px))
    if coeffp[0] < 0 and abs(residp) >= linresid[0]:
        goodP = True
    else:
        goodP = False
        goodS = False
        return sch,residp,resids,spec_px,freq_px,spec_sx,freq_sx,pspec, \
               s1spec,s2spec,pfreq,s1freq,s2freq,pn_spec,s1n_spec, \
               s2n_spec,pn_freq,s1n_freq,s2n_freq,frmin,frmax,goodP,goodS
    
    ## FIND THE LONGEST SEGMENT OF S SPECTRUM WITH SNR > SNRCRT
    if SDATA:
        (begind1,endind1,frmins1,frmaxs1,franges1) = longseg(s1snr,snrcrts,s1freq)
        (begind2,endind2,frmins2,frmaxs2,franges2) = longseg(s2snr,snrcrts,s2freq)
        if franges1 >= franges2:
            begind = begind1
            endind = endind1
            frmins = frmins1
            frmaxs = frmaxs1
            franges = franges1
            sch = 1
        else:
            begind = begind2
            endind = endind2
            frmins = frmins2
            frmaxs = frmaxs2
            franges = franges2
            sch = 2
        frmin[1] = frmins
        frmax[1] = frmaxs
        if franges < snrcrts[1] or frmins > 4:
            goodS = False
        elif sch == 1 :
            goodS = True
            spec_sx = s1spec[begind:endind]
            freq_sx = s1freq[begind:endind]
            spec_nx = s1n_spec[begind:endind]
            nsamp[1] = freq_sx.shape[0]
            snr_sx = s1snr[begind:endind]
            snrmed[1] = float(np.median(snr_sx))
            coeffs = np.polyfit(freq_sx,np.log(spec_sx),1)
            syns = coeffs[1]+freq_sx*coeffs[0]
            resids = seis.lincorrcoef(freq_sx,np.log(spec_sx))
            if coeffs[0] < 0 and abs(resids) >= linresid[1]:
                goodS = True
            else:
                goodS = False
        elif sch == 2 :
            goodS = True
            spec_sx = s2spec[begind:endind]
            freq_sx = s2freq[begind:endind]
            spec_nx = s2n_spec[begind:endind]
            nsamp[1] = freq_sx.shape[0]
            snr_sx = s2snr[begind:endind]
            snrmed[1] = float(np.median(snr_sx))
            coeffs = np.polyfit(freq_sx,np.log(spec_sx),1)
            syns = coeffs[1]+freq_sx*coeffs[0]
            resids = seis.lincorrcoef(freq_sx,np.log(spec_sx))
            if coeffs[0] < 0 and abs(resids) >= linresid[1]:
                goodS = True
            else:
                goodS = False 

    # Plot snr of P and S wave spectra        
    if doplotsnr:
        snr_figure_dir = figdir2+"/%s" %(orid)
        mkdir(snr_figure_dir)
        ## Plot P and S wave
        if SDATA and goodS:
            plt.figure(3)
            plt.clf()
            plt.figure(figsize=(12,6),dpi=200)
            # Plot spectra of P wave signal and noise
            plt.subplot(2,2,1)
            plt.plot(pfreq,np.log(pspec),'b',label='signal')  
            plt.plot(pn_freq,np.log(pn_spec),'r',label='noise')
            plt.plot(freq_px,np.log(spec_px),'k',label='signal in frange')
            plt.plot(freq_px,synp,'g--',linewidth=2,label='synthetic')
            plt.plot([frmin[0],frmin[0]],np.log([min(pn_spec),max(pspec)]),'g')
            plt.plot([frmax[0],frmax[0]],np.log([min(pn_spec),max(pspec)]),'g')
            plt.text(10,max(np.log(pspec))-1,'slope = %.4f' % coeffp[0])
            plt.text(10,max(np.log(pspec))-3,'lincorr = %.4f' % residp)
            plt.xlim([0,50])
            plt.ylabel('ln(Ap) on '+chan[2]+', nm/s')
            plt.title('Station: %s' % sta) 
            plt.legend(loc='upper right')
            # Plot SNR of P wave spectra           
            plt.subplot(2,2,2)
            plt.plot(pfreq,psnr)
            plt.axhline(snrcrtp[0],color='g',linestyle='--')
            plt.axvline(frmin[0],color='g')
            plt.axvline(frmax[0],color='g')
            plt.xlim([0,50])
            plt.ylim([0,1.2*max(psnr)])
            plt.ylabel('P Signal-to-Noise Ratio')
            plt.title('%s' % orid) 
            # Plot spectra of S wave signal and noise
            plt.subplot(2,2,3)
            if sch == 1:
                plt.plot(s1freq,np.log(s1spec),'b',label='signal')  
                plt.plot(s1n_freq,np.log(s1n_spec),'r',label='noise') 
            elif sch == 2:
                plt.plot(s2freq,np.log(s2spec),'b',label='signal')   
                plt.plot(s2n_freq,np.log(s2n_spec),'r',label='noise') 
            plt.plot(freq_sx,np.log(spec_sx),'k',label='signal in frange')
            plt.plot(freq_sx,syns,'g--',linewidth=2,label='synthetic')
            plt.plot([frmin[1],frmin[1]],np.log([min(s1n_spec),max(s1spec)]),'g')
            plt.plot([frmax[1],frmax[1]],np.log([min(s1n_spec),max(s1spec)]),'g')
            plt.xlim([0,50])
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('ln(As) on '+chan[sch-1]+', nm/s')
            plt.text(10,max(np.log(s1spec))-1,'slope = %.4f' % coeffs[0])
            plt.text(10,max(np.log(s1spec))-3,'lincorr = %.4f' % resids)
            plt.legend(loc='upper right')
            # Plot SNR of S wave spectra
            plt.subplot(2,2,4)
            if sch == 1:
                plt.plot(s1freq,s1snr)
                plt.ylim([0,1.2*max(s1snr)])
            elif sch == 2:
                plt.plot(s2freq,s2snr)
                plt.ylim([0,1.2*max(s2snr)])
            plt.xlim([0,50])
            plt.axhline(snrcrts[0],color='g',linestyle='--')
            plt.axvline(frmin[1],color='g')
            plt.axvline(frmax[1],color='g')
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('S Signal-to-Noise Ratio')
            plt.savefig(snr_figure_dir+"/%s_%s_PSsnr.pdf" %(net,sta))
        else:
            ## Plot spectra of only P wave
            plt.figure(4)
            plt.clf()
            # Plot spectra of P wave signal and noise
            plt.subplot(2,1,1)
            plt.plot(pfreq,np.log(pspec),'b',label='signal') 
            plt.plot(pn_freq,np.log(pn_spec),'r',label='noise')
            plt.plot(freq_px,np.log(spec_px),'k',label='signal in frange')
            plt.plot(freq_px,synp,'g--',linewidth=2,label='synthetic')
            plt.plot([frmin[0],frmin[0]],np.log([min(pn_spec),max(pspec)]),'g')
            plt.plot([frmax[0],frmax[0]],np.log([min(pn_spec),max(pspec)]),'g')
            plt.text(10,max(np.log(pspec))-1,'slope = %.4f' % coeffp[0])
            plt.text(10,max(np.log(pspec))-3,'lincorr = %.4f' % residp)
            plt.xlim([0,50])
            plt.ylabel('ln(Ap) on '+chan[2]+', nm/s')
            plt.title('Station: %s' % sta) 
            plt.legend(loc='upper right')
            # Plot SNR of P wave spectra           
            plt.subplot(2,1,2)
            plt.plot(pfreq,psnr)
            plt.axhline(snrcrtp[0],color='g',linestyle='--')
            plt.axvline(frmin[0],color='g')
            plt.axvline(frmax[0],color='g')
            plt.xlim([0,50])
            plt.ylim([0,1.2*max(psnr)])
            plt.ylabel('P Signal-to-Noise Ratio')
            plt.title('%s' % orid)
            plt.savefig(snr_figure_dir+"/%s_%s_Psnr.pdf" %(net,sta))
            
    return sch,residp,resids,spec_px,freq_px,spec_sx,freq_sx,pspec, \
           s1spec,s2spec,pfreq,s1freq,s2freq,pn_spec,s1n_spec, \
           s2n_spec,pn_freq,s1n_freq,s2n_freq,frmin,frmax,goodP,goodS


def buildG(saving,stalst,alpha,POS,icase,mw_flag):
##  Build G matrix
##          G = [[1, -pi*f1i*f1i**(-alpha), 0, ..., 0],         ##
##               [1, 0, -pi*f2i*f2i**(-alpha), ..., 0],         ##
##               [1, 0, 0, ..., -pi*fMi*fMi**(-alpha)]]         ##
##  INPUT:  (saving) ==> saved frequency for each station: saving[sta][1]['p']
##          (stalst) ==> list of used stations
##          (alpha)  ==> alpha value(s)
##          (POS)    ==> 'P' or 'S'
##          (icase)  ==> 1: high quality for finding best fc and alpha
##                       2: low quality for t* inversion
##                       3: low quality for t* inversion without bad fitting
##  OUTPUT: (G)      ==> G matrix for t* inversion
    if icase == 3:
        icase = 2
    if POS.upper() == 'P':
        ind = 0
    elif POS.upper() == 'S':
        ind = 1
    else:
        raise ValueError("P or S wave?")
    for ista in range(len(stalst)):
        sta = stalst[ista]
        for alco in range(len(alpha)):
            freq_x = saving[sta][icase][POS.lower()][0]
            exponent = -1*np.pi*freq_x*(freq_x**(-alpha[alco]))
            exponent = np.array([exponent]).transpose()
            if alco == 0:
                Gblock = np.atleast_3d(exponent)
            else:
                Gblock = np.dstack((Gblock,exponent))
        if ista == 0:
            G = Gblock
        else:
            oldblock = np.hstack((G,np.zeros((G.shape[0],1,len(alpha)))))
            newblock = np.hstack((np.zeros((Gblock.shape[0],G.shape[1],len(alpha))),Gblock))
            G = np.vstack((oldblock,newblock))
    if mw_flag == False:
        if POS.upper() == 'P':
            G = np.hstack((np.ones((G.shape[0],1,len(alpha))),G))
    
    return G


def buildd(saving,stalst,fc,POS,icase,lnM=0):
##  Build data matrix
##          d = [ln(A1)-ln(C1)+ln(1+(f1i/fc)**2),                   ##
##               ln(A2)-ln(C2)+ln(1+(f2i/fc)**2),                   ##
##               ln(AM)-ln(CM)+ln(1+(fMi/fc)**2)]                   ##
##  INPUT:  (saving) ==> saved spectrum for each station: saving[sta][1]['p']
##          (stalst) ==> list of used stations
##          (fc)     ==> corner frequency
##          (POS)    ==> 'P' or 'S'
##          (icase)  ==> 1: high quality for finding best fc and alpha
##                       2: low quality for t* inversion
##                       3: low quality for t* inversion without bad fitting
##          (lnM)    ==> when POS='S', log of seismic moment
##  OUTPUT: (data)   ==> data matrix for t* inversion
    if icase == 3:
        icase = 2
    if POS.upper() == 'P':
        ind = 0
    elif POS.upper() == 'S':
        ind = 1
    else:
        raise ValueError("P or S wave?")
    for ista in range(len(stalst)):
        sta = stalst[ista]
        freq_x = saving[sta][icase][POS.lower()][0]
        spec_x = saving[sta][icase][POS.lower()][1]
        correc = saving[sta]['corr'][ind]
        stad = np.array([np.log(spec_x)-np.log(correc)
                        +np.log(1+(freq_x/fc)**2)-lnM]).transpose()
        if ista == 0:
            data = stad
        else:
            data = np.vstack((data,stad))
    
    return data


def calresspec(saving,sta,orid,POS,lnmomen,fc,alpha):
## CALCULATE RESIDUAL SPECTRUM FOR EACH STATION
    if POS.upper() == 'P':
        ind = 0
    elif POS.upper() == 'S':
        ind = 1
    else:
        raise ValueError("P or S wave?")
    freq_x = saving[2][POS.lower()][0]
    
    spec_x = saving[2][POS.lower()][1]
    correc = saving['corr'][ind]
    invtstar = saving[3]['tstar'][ind]
    righthand = lnmomen-np.pi*freq_x*(freq_x**(-alpha)*invtstar)
    resspec = np.array([np.log(spec_x)-np.log(correc)
                  +np.log(1+(freq_x/fc)**2)-righthand])
    resratio = resspec/righthand*100
    resspec = np.vstack((freq_x,resspec))
    resspec = np.vstack((resspec,resratio))
    resspec = resspec.transpose()
    return resspec


def fitting(saving,sta,orid,POS,lnmomen,fc,alpha,icase):
## CALCULATE HOW WELL THE SYNTHETIC SPECTRUM FITS THE DATA
## IF THE FITTING CURVE IS BELOW THE NOISE, THEN resid = 999999.
    if POS.upper() == 'P':
        ind = 0
        sch = 0
    elif POS.upper() == 'S':
        ind = 1
        sch = saving['sch']
    else:
        raise ValueError("P or S wave?")
    corr = saving['corr'][ind]
    spec = saving['spec'][sch]
    freq = saving['freq'][sch]
    frmin = saving[2]['frmin'][ind]
    frmax = saving[2]['frmax'][ind]
    invtstar = saving[icase]['tstar'][ind]
    synspec = (corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar)/(1+(freq/fc)**2))
    indx = np.all([(freq>=frmin),(freq<frmax)],axis=0)
    specx = spec[indx]
    freqx = freq[indx]
    synx = synspec[indx]
    resid=(1-((np.linalg.norm(np.log(synx)-np.log(specx)))**2/(len(freqx)-1)
             /np.var(np.log(specx))))    
    ##    RESID = 1-√(∑((ln(A(synthetic))-ln(A(observed)))^2))/(n-1)/σ(ln(A(observed))) Bill Menke, 'geophysical data analysis discrete inverse theory'
    resid1 = (np.linalg.norm(np.log(synx)-np.log(specx))/(len(freqx)))/0.30
    ##    L2-NORM MISFIT = √(∑((ln(A(synthetic))-ln(A(observed)))^2))/n. Here 0.30 is just for normalization
    resid2 = (np.linalg.norm(np.log(synx)-np.log(specx),ord=1)/(len(freqx)))/1.50
    ##    L1-NORM MISFIT = (∑|ln(A(synthetic))-ln(A(observed))|)/n. Here 1.50 is just for normalization
    resid3 = (1-2*np.sum(np.log(synx)*np.log(specx))/(np.sum((np.log(specx))**2)+np.sum((np.log(synx))**2)))/0.8
    ##    CORRELATIVE FUNCTION. 1-2*∑(ln(A(synthetic))*ln(A(observed)))/∑((ln(A(synthetic))^2)+(ln(A(observed)))^2). Here 0.80 is just for normalization
    df = abs(freq[1]-freq[0])
    nlowf = 0
    narea = 0

    return resid,resid1,resid2,resid3


def plotspec(plotspecloglog,fitting,fitting1,fitting2,fitting3,lincorr,saving,net,sta,orid,chan,POS,lnmomen,fc,alpha,icase,figdir5,sitedata=0):
## PLOT AMPLITUDE SPECTRUM
    if POS.upper() == 'P':
        ind = 0
        sch = 0
        xmax = 20
        ymin = -10
        ymax = 10
        textx = 6
    elif POS.upper() == 'S':
        ind = 1
        sch = saving['sch']
        xmax = 20
        ymin = -10
        ymax = 10
        textx = 2.5
    else:
        raise ValueError("P or S wave?")
    plot_chan = [chan[2], chan[0],chan[1]]
    corr = saving['corr'][ind]
    spec = saving['spec'][sch]
    freq = saving['freq'][sch]
    n_spec = saving['nspec'][ind]
    n_freq = saving['nfreq'][ind]
    frmin = saving[2]['frmin'][ind]
    frmax = saving[2]['frmax'][ind]
    invtstar = saving[icase]['tstar'][ind]
    synspec = (corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar)/(1+(freq/fc)**2))
    if POS.upper() == 'S':
        invtstarP = saving[icase]['tstar'][0]
        ttP = saving['Ptt']
        ttS = saving['Stt']
        QpQs = 2.25
        invtstar2 = invtstarP*QpQs*ttS/ttP
        synspec2 = (corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar2)/(1+(freq/fc)**2))
        QpQs = 1.75
        invtstar2 = invtstarP*QpQs*ttS/ttP
        synspec3 = (corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar2)/(1+(freq/fc)**2))
    indx = np.all([(freq>=frmin),(freq<frmax)],axis=0)
    specx = spec[indx]
    freqx = freq[indx]
    synx = synspec[indx]
    resid = (1-((np.linalg.norm(np.log(synx)-np.log(specx)))**2/(len(freqx)-1)
             /np.var(np.log(specx))))    
    df = abs(freq[1]-freq[0])
    nlowf = 0
    narea = 0
    for ifreq in range(len(freq)):
        if (freq[ifreq]>frmax and freq[ifreq]<15):
            if (np.log(synspec[ifreq])<np.log(spec[ifreq]) or
                 np.log(synspec[ifreq])>np.log(spec[ifreq])+1):
                narea = narea+np.log(spec[ifreq])-np.log(synspec[ifreq])
            if np.log(synspec[ifreq])>np.log(spec[ifreq])+2:
                nlowf = nlowf+5
            elif np.log(synspec[ifreq])>np.log(spec[ifreq])+1:
                nlowf = nlowf+1
    if narea<-10 and nlowf*df>3:
        resid = 0
    
    #####################remain to be studied########################
    if not isinstance(sitedata, int):
        plt.plot(freq,np.log(spec)-sitedata,'b--')   
    #################################################################
    if plotspecloglog:
        spec_figure_dir = figdir5+"/%s" %(orid)
        mkdir(spec_figure_dir)
        # X-log axis, Y-log axis. 
        fig = plt.figure(7)
        ax = fig.add_subplot(111)
        plt.clf()
        plt.loglog(freq,spec,'b',label='signal')
        plt.loglog(n_freq,n_spec,'r',label='noise')
        plt.loglog([frmin,frmin],[min(n_spec),max(spec)],'g',label='frmin')
        plt.loglog([frmax,frmax],[min(n_spec),max(spec)],'g',label='frmax')
        plt.loglog(freq,synspec,'g--',linewidth=2,label='synthetic')
        plt.legend(loc='lower left')
        text = "t* = %.2f\n" %(invtstar) + \
               "fitting1 = %.4f\n" %(fitting)+ \
               "fitting2 = %.4f\n" %(fitting1)+ \
               "fitting3 = %.4f\n" %(fitting2)+ \
               "fitting4 = %.4f\n" %(fitting3)+ \
               "frange = [%.4f, %.4f]\n" %(frmin, frmax)+ \
               "sd = %.4f" %(saving[icase]['err'][ind])
        
        plt.text(0.75,0.6,text,transform=ax.transAxes)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('A%s on %s, nm/s' %(POS.lower(),plot_chan[sch]))
        plt.title('%s Station: %s' % (orid,sta))
        plt.savefig(spec_figure_dir+'/%s_%s_%sspectrum_loglog.pdf' % (net,sta,POS.upper()))
        
        return
    else:
        spec_figure_dir = figdir5+"/%s" %(orid)
        mkdir(spec_figure_dir)
        # X-linear axis, Y-log axis.
        fig = plt.figure(8)
        ax = fig.add_subplot(111)
        plt.clf()
        plt.plot(freq,np.log(spec),'b',label='signal')
        plt.plot(n_freq,np.log(n_spec),'r',label='noise')
        plt.plot([frmin,frmin],np.log([min(n_spec),max(spec)]),'g',label='frmin')
        plt.plot([frmax,frmax],np.log([min(n_spec),max(spec)]),'g',label='frmax')
        plt.plot(freq,np.log(synspec),'g--',linewidth=2,label='synthetic')
        plt.legend(loc='lower left')
        text = "t* = %.2f\n" %(invtstar) + \
               "fitting1 = %.4f\n" %(fitting)+ \
               "fitting2 = %.4f\n" %(fitting1)+ \
               "fitting3 = %.4f\n" %(fitting2)+ \
               "fitting4 = %.4f\n" %(fitting3)+ \
               "frange = [%.4f, %.4f]\n" %(frmin, frmax)+ \
               "sd = %.4f" %(saving[icase]['err'][ind])

        plt.text(0.75,0.6,text,transform=ax.transAxes)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('ln(A%s) on %s, nm/s' %(POS.lower(),plot_chan[sch]))
        plt.title('%s  Station: %s.%s' % (orid,net,sta))
        plt.savefig(spec_figure_dir+'/%s_%s_%sspectrum_loglinear.pdf' % (net,sta,POS.upper()))
        
        return