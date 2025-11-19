#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Sub functions for t* inversion 


import os
import glob
import numpy as np
import seissign as seis
import matplotlib as mpl
import multitaper as mtm
from scipy import signal
import tstar_parameters as tp
import matplotlib.pyplot as plt


# Output Type 42 (Truetype) or Type 3 (Open type)
mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['font.size'] = 10
# Use scientific notation if log10 of the axis range is smaller than the first or larger than the second
mpl.rcParams['axes.formatter.limits'] = [-2, 2]



def readseismo(orid, ARRIV):
    """
    Read seismograms in txt format from './data/correctedSeismograms/'.
    """

    dt = ARRIV['dt']
    net = ARRIV['net']
    sta = ARRIV['sta']
    chan = ARRIV['chan']
    pretime = ARRIV['T0']

    for ichan in range(len(chan)):
        sacfile = glob.glob("%s/%s/%s.%s*%s.txt" %(tp.sacdir, orid, net, sta, chan[ichan]))[0]
        if not os.path.isfile(sacfile):
            print('ERROR: %s does not exist' %(sacfile))
            dd = np.array(range(18000))
            tt = np.array(range(18000))
            flag = False
            return dd, tt, flag
        else:
            ddchan = np.fromstring(("".join(open(sacfile).readlines()[30:])), sep=' ')
            nn = ddchan.size
            if ichan == 0:
                nmin = nn
            else:
                if nn < nmin:
                    nmin = nn

    for ichan in range(len(chan)):
        sacfile = glob.glob("%s/%s/%s.%s*%s.txt" %(tp.sacdir, orid, net, sta, chan[ichan]))[0]
        ddchan = np.fromstring(("".join(open(sacfile).readlines()[30:])), sep=' ')
        ddchan = ddchan[:nmin]
        if ichan == 0:
            dd = ddchan
        else:
            dd = np.vstack((dd, ddchan))
    
    flag = True
    tt = -pretime+dt*np.array(range(nmin)) 
    
    return dd, tt, flag


def fixwin(dd, tt, param, ARRIV, orid):
    """
    Window seismic data in time domain.

    """

    dt = ARRIV['dt']
    samplerate = int(1/float(dt))
    net = ARRIV['net']
    sta = ARRIV['sta']
    chan = ARRIV['chan']
    pchan = ARRIV['pchan']

    WLP = param['WLP']
    WLS = param['WLS']
    prewin = param['prewin']
    gaps = param['gaps']
    
    ## Window P wave
    pch = chan.index(pchan)
    pind = np.all([(tt >= -1*prewin[0]), (tt <= (WLP-prewin[0]))], axis=0)
    p_dd = signal.detrend(dd[pch][pind]-np.mean(dd[pch][pind])) # type: ignore
    pnind = np.all([(tt <= -1*prewin[0]), (tt >= (-WLP-prewin[0]))], axis=0)
    pn_dd = signal.detrend(dd[pch][pnind]-np.mean(dd[pch][pnind])) # type: ignore
    
    ## Make sure P wave and noise have same data points
    size = min(p_dd.size, pn_dd.size)
    p_dd = p_dd[0:size]
    pn_dd = pn_dd[0:size]

    ## Plot P wave seismogram
    if param['doplotseis']:
        seis_figure_dir = tp.figdir1+"/%s" %(orid)
        if not os.path.exists(seis_figure_dir):
            os.makedirs(seis_figure_dir)
        plt.figure(1, figsize=(10, 5))
        plt.clf()
        b, a = signal.butter(4, 1, btype="highpass", fs=samplerate)
        filt_dd = signal.filtfilt(b, a, dd[pch])
        plt.plot(tt, filt_dd)
        plt.xlim([-5, 5])
        P_begind = int((-5-tt[0])/dt)
        P_endind = int((5-tt[0])/dt)
        plot_dd = filt_dd[P_begind:P_endind]
        plt.ylim(-1.2*max(abs(plot_dd)), 1.2*max(abs(plot_dd)))
        plt.xlabel("Time relative to P arrival time (s)")
        plt.ylabel("Velocity Amplitude (nm/s)")
        plt.axvline(-1*prewin[0], color="green", linewidth=1.0)
        plt.axvline(WLP-prewin[0], color="green", linewidth=1.0)
        plt.axvline(-WLP-prewin[0], color="green", linewidth=1.0)
        plt.text(-WLP-prewin[0], -1.1*max(abs(plot_dd)), "Noise", fontsize=10)
        plt.text(-1*prewin[0], -1.1*max(abs(plot_dd)), "P wave", fontsize=10)
        figure_title = "%s %s.%s.%s" %(orid, net, sta, pchan)
        plt.title(figure_title)
        plt.savefig(seis_figure_dir+'/%s_%s_Pwindata.pdf' % (net, sta))
    
    ## Window S wave on both horizontal channels
    if ARRIV['SDATA']:
        sminusp = ARRIV['T1'] - ARRIV['T0']
        if sminusp < (WLS+WLP+prewin[0]-prewin[1]+gaps):
            print('P & S arrivels are too close - proceeding as if no S pick')
            ARRIV['SDATA'] = False
        if round(sminusp/dt,5) == int(sminusp/dt):
            sminusp += 0.001

        ## Horizontal channel 1
        sch = 0
        sind = np.all([(tt>=(sminusp-prewin[1])),(tt<=(sminusp+WLS-prewin[1]))],axis=0)
        s_dd1 = signal.detrend(dd[sch][sind]-np.mean(dd[sch][sind]))
        ## Noise is defined as gaps s before S arrival
        snind = np.all([(tt<=(sminusp-prewin[1]-gaps)),(tt>=(sminusp-WLS-prewin[1]-gaps))],axis=0)
        sn_dd1 = signal.detrend(dd[sch][snind]-np.mean(dd[sch][snind]))
        ## P coda is defined as right before S arrival
        pcind = np.all([(tt<=(sminusp-prewin[1])),(tt>=(sminusp-WLS-prewin[1]))],axis=0)
        pc_dd1 = signal.detrend(dd[sch][pcind]-np.mean(dd[sch][pcind]))
        ## Make sure S wave, noise and P code wave have same data points
        minlen = min(s_dd1.size, sn_dd1.size, pc_dd1.size)
        s_dd1 = s_dd1[0:minlen]
        sn_dd1 = sn_dd1[0:minlen]
        pc_dd1 = pc_dd1[0:minlen]

        ## Horizontal channel 2
        sch = 1
        sind = np.all([(tt>=(sminusp-prewin[1])),(tt<=(sminusp+WLS-prewin[1]))],axis=0)
        s_dd2 = signal.detrend(dd[sch][sind]-np.mean(dd[sch][sind]))
        ## Noise is defined as gaps s before S arrival
        snind = np.all([(tt<=(sminusp-prewin[1]-gaps)),(tt>=(sminusp-WLS-prewin[1]-gaps))],axis=0)
        sn_dd2 = signal.detrend(dd[sch][snind]-np.mean(dd[sch][snind]))
        ## P coda is defined as right before S arrival
        pcind = np.all([(tt<=(sminusp-prewin[1])),(tt>=(sminusp-WLS-prewin[1]))],axis=0)
        pc_dd2 = signal.detrend(dd[sch][pcind]-np.mean(dd[sch][pcind]))
        ## Make sure S wave, noise and P code wave have same size
        minlen = min(s_dd2.size, sn_dd2.size, pc_dd2.size)
        s_dd2 = s_dd2[0:minlen]
        sn_dd2 = sn_dd2[0:minlen]
        pc_dd2 = pc_dd2[0:minlen]

        ## Plot S wave seismogram
        if param['doplotseis']:
            ## Plot seismogram for horizontal channel 1
            plt.figure(2, figsize=(10, 10))
            plt.clf()
            plt.subplot(2, 1, 1)
            b, a = signal.butter(4, 1, btype="highpass", fs=samplerate)
            filt_dd1 = signal.filtfilt(b, a, dd[0])
            plt.plot(tt, filt_dd1)
            plt.xlim([sminusp-10, sminusp+10])
            S_begind = int((sminusp-10-tt[0])/dt)
            S_endind = int((sminusp+10-tt[0])/dt)
            plot_S_dd1 = filt_dd1[S_begind:S_endind]
            y_max = max(abs(plot_S_dd1))
            plt.ylim(-1.2*y_max, 1.2*y_max)
            plt.axvline(sminusp-prewin[1], color="green")
            plt.axvline(sminusp+WLS-prewin[1], color="green")
            plt.axvline(sminusp-prewin[1]-gaps, color="green")
            plt.axvline(sminusp-WLS-prewin[1]-gaps, color="green")
            plt.axvline(-1*prewin[1], color="green")
            plt.text(sminusp-WLS-prewin[1]-gaps, -1.1*y_max, "Noise", fontsize=10)
            plt.text(sminusp-prewin[1], -1.1*y_max, "S wave", fontsize=10)
            plt.ylabel("Velocity Amplitude (nm/s)")
            figure_title = "%s %s.%s.%s" %(orid, net, sta, chan[0])
            plt.title(figure_title)
            
            ## Plot seismogram for horizontal channel 2
            plt.subplot(2, 1, 2)
            filt_dd2 = signal.filtfilt(b, a, dd[1])
            plt.plot(tt, filt_dd2)
            plt.xlim([sminusp-10,sminusp+10])
            S_begind = int((sminusp-10-tt[0])/dt)
            S_endind = int((sminusp+10-tt[0])/dt)
            plot_S_dd2 = filt_dd2[S_begind:S_endind]
            y_max = max(abs(plot_S_dd2))
            plt.ylim(-1.2*y_max, 1.2*y_max)
            plt.axvline(sminusp-prewin[1], color="green")
            plt.axvline(sminusp+WLS-prewin[1], color="green")
            plt.axvline(sminusp-prewin[1]-gaps, color="green")
            plt.axvline(sminusp-WLS-prewin[1]-gaps, color="green")
            plt.axvline(-1*prewin[1], color="green")
            plt.text(sminusp-WLS-prewin[1]-gaps, -1.1*y_max, "Noise", fontsize=10)
            plt.text(sminusp-prewin[1], -1.1*y_max, "S wave", fontsize=10)
            plt.xlabel("Time relative to P arrival time (s)")
            plt.ylabel("Velocity Amplitude (nm/s)")
            figure_title = "%s %s.%s.%s" %(orid, net, sta, chan[1])
            plt.title(figure_title)
            plt.savefig(seis_figure_dir+'/%s_%s_Swindata.pdf' % (net, sta))
    
    else:
        s_dd1 = p_dd
        sn_dd1 = pn_dd
        pc_dd1 = p_dd
        s_dd2 = p_dd
        sn_dd2 = pn_dd
        pc_dd2 = p_dd

    return p_dd, pn_dd, s_dd1, sn_dd1, pc_dd1, s_dd2, sn_dd2, pc_dd2


def longseg(dt, snr, snrcrtpara, freq, param):
    """
    ## Find the longest segment of spectra with SNR >= threshold
    """

    samplerate = int(1/float(dt))
    for samples in [40, 50, 100]:
        if samplerate == samples:
            minf = param[str(samples)][0]
            maxf = param[str(samples)][1]
    
    lenspec = len([ifreq for ifreq in freq if (ifreq<maxf and ifreq>=minf)])
    ind1 = int(min(np.nonzero(freq>=minf)[0]))
    ind2 = int(max(np.nonzero(freq<maxf)[0]))
    w = 0
    m = []
    bindex = []
    eindex = []
    snrcrt = snrcrtpara[0]
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
        
        return begind, endind, frmin, frmax, frange
    

    ## Find the longest segment
    longest = m.index(max(m))
    frmin = freq[bindex[longest]]
    frmax = freq[eindex[longest]]
    frange = frmax-frmin
    
    ## Favor the second longest segment if it has lower frequency < 4 Hz and longer than 1/4 of the longest one
    if len(m) >= 2:
        for mind in list(reversed(range(len(m)))):
            mii = mind-len(m)
            longest2 = m.index(sorted(m)[mii])
            frmin2 = freq[bindex[longest2]]
            if frmin2 <= 2.0:
                frmax2 = freq[eindex[longest2]]
                frange2 = frmax2-frmin2
                if frmin2<frmin and 4*frange2>frange and frange2>snrcrtpara[1]:
                    frmin = frmin2
                    frmax = frmax2
                    frange = frange2
                    longest = longest2
                    break
    
    begind = bindex[longest]
    endind = eindex[longest]
    
    ## Extend frequency band to low SNR
    if snrcrtpara[2] < snrcrtpara[0]:
        if begind > ind1+1:
            while snr[begind-1]<snr[begind] and snr[begind-1]>snrcrtpara[2] and begind-1>ind1+1:
                begind = begind-1
        if endind < ind2-1:
            while snr[endind+1]<snr[endind] and snr[endind+1]>snrcrtpara[2] and endind+1<ind2-1:
                endind = endind+1
        frmin = freq[begind]
        frmax = freq[endind]
        frange = frmax-frmin
    
    return begind,endind,frmin,frmax,frange


def dospec(pwindata, swindata1, swindata2, orid, param, data_quality, ARRIV, ORIG):
    """
    Calculate amplitude spectrum of windowed waveform using multi-taper method
    """
    
    dt = ARRIV['dt']
    net = ARRIV['net']
    sta = ARRIV['sta']
    chan = ARRIV['chan']

    snrcrtpara_p = param['snrcrtp'+str(data_quality)]
    snrcrtpara_s = param['snrcrts'+str(data_quality)]
    smlen = 11
    spec = []
    freq = []
    n_spec = []
    n_freq = []

    ## Determine P wave spectra
    for ii in range(pwindata.shape[0]):
        if param['mtspec_para'] == 1:
            nft = 1024
            mtmresult = mtm.mtspec(pwindata[ii], dt, time_bandwidth=5, nfft=nft)
        elif param['mtspec_para'] == 2:
            mtmresult = mtm.sine_psd(pwindata[ii], dt)
        
        newspec = mtmresult[0][1:]
        newfreq = mtmresult[1][1:]
        ## Convert velocity to displacement by dividing by 2*pi*f (Gubbins, p30)
        newspec = np.sqrt(newspec)/(2*np.pi*newfreq)
        if ii == 0:
            spec.append(newspec)
            freq.append(newfreq)
        else:
            n_spec.append(newspec)
            n_freq.append(newfreq)
 
    ## Determine S wave spectra on channel 1
    for ii in range(swindata1.shape[0]):
        if param['mtspec_para'] == 1:
            nft = 1024
            mtmresult = mtm.mtspec(swindata1[ii], dt, time_bandwidth=5, nfft=nft)
        elif param['mtspec_para'] == 2:
            mtmresult = mtm.sine_psd(swindata1[ii], dt)
        
        newspec = mtmresult[0][1:]
        newfreq = mtmresult[1][1:]
        newspec = np.sqrt(newspec)/(2*np.pi*newfreq)

        if ii == 0:
            spec.append(newspec)
            freq.append(newfreq)
        elif ii == 1:
            n_spec.append(newspec)
            n_freq.append(newfreq)
        elif ii == 2:
            pcspec = newspec
    
    ## Determine S wave spectra on channel 2
    for ii in range(swindata2.shape[0]):
        if param['mtspec_para'] == 1:
            nft = 1024
            npi = 3.0 # not used in following codes
            mtmresult = mtm.mtspec(swindata2[ii], dt, time_bandwidth=5, nfft=nft)
        elif param['mtspec_para'] == 2:
            mtmresult = mtm.sine_psd(swindata2[ii], dt)
        
        newspec = mtmresult[0][1:]
        newfreq = mtmresult[1][1:]
        newspec = np.sqrt(newspec)/(2*np.pi*newfreq)

        if ii == 0:
            spec.append(newspec)
            freq.append(newfreq)
        elif ii == 1:
            n_spec.append(newspec)
            n_freq.append(newfreq)
    
    spec_px = spec[0]
    freq_px = freq[0]
    spec_sx = spec[1]
    freq_sx = freq[1]
    frmin = [6, 6]
    frmax = [6, 6]
    goodP = False
    goodS = False
    
    # SNR of P wave spectra
    p_spec = spec[0]
    p_freq = freq[0]
    pn_spec = n_spec[0]
    pn_freq = n_freq[0]
    psnr = p_spec/pn_spec
    psnr = seis.smooth(p_spec, smlen)/seis.smooth(pn_spec, smlen)
    
    (begindp, endindp, frminp, frmaxp, frangep) = longseg(dt, psnr, snrcrtpara_p, p_freq, param)
    frmin[0] = frminp
    frmax[0] = frmaxp
    freq_px = p_freq[begindp:endindp]
    spec_px = p_spec[begindp:endindp]
    
    try:
        coeffp = np.polyfit(freq_px, np.log(spec_px), 1)
        synp = coeffp[1]+freq_px*coeffp[0]
    except:
        coeffp = [0.1, 0.1]
    
    residp = seis.lincorrcoef(freq_px, np.log(spec_px))

    if param['doplotsnr']:
        # Plot P wave spectra
        snr_figure_dir = tp.figdir2+"/%s" %(orid)
        if not os.path.exists(snr_figure_dir):
            os.makedirs(snr_figure_dir)
        plt.figure(3, figsize=(8, 8))
        plt.clf()
        # Plot SNR of P wave spectra
        plt.subplot(2, 1, 1)
        plt.plot(p_freq, psnr)
        plt.axhline(snrcrtpara_p[0], color='g', linestyle='--')
        plt.axvline(frmin[0], color='g')
        plt.axvline(frmax[0], color='g')
        plt.xlim([0, int(max(p_freq)+1)])
        plt.ylim([0, 1.2*max(psnr)])
        plt.ylabel('P Signal-to-Noise Ratio')
        plt.title('%s, lat: %s, lon: %s, depth: %s' %(orid, ORIG['lat'], ORIG['lon'], ORIG['dep']))
        # Plot spectra of P wave signal and noise
        plt.subplot(2, 1, 2)
        plt.plot(p_freq, np.log(p_spec),'b',label='signal') 
        plt.plot(pn_freq, np.log(pn_spec),'r',label='noise')
        plt.plot([frmin[0], frmin[0]], np.log([min(pn_spec), max(p_spec)]),'g')
        plt.plot([frmax[0], frmax[0]], np.log([min(pn_spec), max(p_spec)]),'g')
        plt.text(10, max(np.log(p_spec))-1, 'frmin = %.4f' % frminp)
        plt.text(10, max(np.log(p_spec))-2, 'frmax = %.4f' % frmaxp)
        plt.xlim([0, int(max(p_freq)+1)])
        plt.ylabel('ln(Ap) on '+chan[2]+', nm/s')
        sampling_rate = int(1/dt)
        plt.title('Station: %s, Sampling rate: %s Hz' % (sta, sampling_rate))
        plt.plot(freq_px, np.log(spec_px), 'k', label='signal in frange')
        try:
            plt.plot(freq_px, synp, 'g--', linewidth=2, label='synthetic')
            plt.text(10, max(np.log(p_spec))-3, 'slope = %.4f' % coeffp[0])
            plt.text(10, max(np.log(p_spec))-4, 'residual = %.4f' % residp)
        except:
            pass
        plt.legend(loc='upper right')           
        plt.savefig(snr_figure_dir+"/%s_%s_Psnr.pdf" %(net,sta))

    if frangep < snrcrtpara_p[1] or frminp > param['frminp']:
        goodP = False
        goodS = False
        return goodP, goodS, spec, freq, n_spec, n_freq, spec_px, freq_px, spec_sx, freq_sx, frmin, frmax
    else:
        goodP = True
    
    if coeffp[0] < 0 and abs(residp) >= param['lincor'][0]:
        goodP = True
    else:
        goodP = False
        goodS = False
        
        return goodP, goodS, spec, freq, n_spec, n_freq, spec_px, freq_px, spec_sx, freq_sx, frmin, frmax
    

    if ARRIV['SDATA']:
        # Channel 1
        s1_spec = spec[1]
        s1_freq = freq[1]
        s1n_spec = n_spec[1]
        s1n_freq = n_freq[1]
        s1snr = s1_spec/s1n_spec
        s1snr = seis.smooth(s1_spec,smlen)/seis.smooth(s1n_spec,smlen)
        (beginds1, endinds1, frmins1, frmaxs1, franges1) = longseg(dt, s1snr, snrcrtpara_s, s1_freq, param)
        
        # Channel 2
        s2_spec = spec[2]
        s2_freq = freq[2]
        s2n_spec = n_spec[2]
        s2n_freq = n_freq[2]
        s2snr = s2_spec/s2n_spec
        s2snr = seis.smooth(s2_spec,smlen)/seis.smooth(s2n_spec,smlen)
        (beginds2, endinds2, frmins2, frmaxs2, franges2) = longseg(dt, s2snr, snrcrtpara_s, s2_freq, param)

        if franges2 >= franges1:
            frmins = frmins2
            frmaxs = frmaxs2
            beginds = beginds2
            endinds = endinds2
            franges = franges2
            frmin[1] = frmins
            frmax[1] = frmaxs
            freq_sx = s2_freq[beginds:endinds]
            spec_sx = s2_spec[beginds:endinds]
        else:
            frmins = frmins1
            frmaxs = frmaxs1
            beginds = beginds1
            endinds = endinds1
            franges = franges1
            frmin[1] = frmins
            frmax[1] = frmaxs
            freq_sx = s1_freq[beginds:endinds]
            spec_sx = s1_spec[beginds:endinds]
        
        try:
            coeffs = np.polyfit(freq_sx, np.log(spec_sx), 1)
            syns = coeffs[1]+freq_sx*coeffs[0]
        except:
            coeffs = [0.1, 0.1]
        
        residus = seis.lincorrcoef(freq_sx, np.log(spec_sx))

        if param['doplotsnr']:
            plt.figure(4, figsize=(16, 16))
            plt.clf()
            # Channel 1
            plt.subplot(2, 2, 1)
            plt.plot(s1_freq, s1snr)
            plt.axhline(snrcrtpara_s[0], color='g', linestyle="--")
            plt.axhline(snrcrtpara_s[2], color='g', linestyle="--")
            plt.axvline(frmins, color='g')
            plt.axvline(frmaxs, color='g')
            plt.xlim([0, int(max(s1_freq))+1])
            plt.ylabel("S Signal-to-noise Ratio Channel 1")
            plt.title('%s, lat: %s, lon: %s, depth: %s' %(orid,ORIG['lat'], ORIG['lon'], ORIG['dep']))
            plt.subplot(2, 2, 2)
            plt.plot(s1_freq, np.log(s1_spec), 'b', label='signal')
            plt.plot(s1n_freq, np.log(s1n_spec), 'r', label='noise')
            plt.plot([frmins,frmins], np.log([min(s1n_spec), max(s1_spec)]),'g')
            plt.plot([frmaxs,frmaxs], np.log([min(s1n_spec), max(s1_spec)]),'g')
            plt.text(10, max(np.log(s1_spec))-1, 'frmin = %.4f' % frmins)
            plt.text(10, max(np.log(s1_spec))-2, 'frmax = %.4f' % frmaxs)
            plt.xlim([0, int(max(s1_freq)+1)])
            plt.ylabel("ln(As) on" +chan[0]+' nm/s')
            sampling_rate = int(1/dt)
            plt.title('Station: %s, Sampleing rate: %s Hz' %(sta, sampling_rate))
            if franges1 > franges2:
                plt.text(15, max(np.log(s1_spec))-1, 'prefered channel')
            
            plt.plot(freq_sx, np.log(spec_sx), 'k', label='signal in frange')
            try:
                plt.plot(freq_sx,syns, 'g--', linewidth=2, label='synthetic')
                plt.text(10, max(np.log(s1_spec))-3, 'slope = %.4f' % coeffs[0])
                plt.text(10, max(np.log(s1_spec))-4, 'residual = %.4f' % residus)
            except:
                pass
            plt.legend(loc='upper right')

            # Channel 2
            plt.subplot(2, 2, 3)
            plt.plot(s2_freq, s2snr)
            plt.axhline(snrcrtpara_s[0], color='g', linestyle="--")
            plt.axhline(snrcrtpara_s[2], color='g', linestyle="--")
            plt.axvline(frmins, color='g')
            plt.axvline(frmaxs, color='g')
            plt.xlim([0, int(max(s2_freq))+1])
            plt.ylabel("S Signal-to-noise Ratio Channel 2")
            plt.title('%s, lat: %s, lon: %s, depth: %s' %(orid,ORIG['lat'],ORIG['lon'],ORIG['dep']))
            plt.subplot(2, 2, 4)
            plt.plot(s2_freq, np.log(s2_spec), 'b', label='signal')
            plt.plot(s2n_freq, np.log(s2n_spec), 'r', label='noise')
            plt.plot([frmins,frmins], np.log([min(s2n_spec), max(s2_spec)]), 'g')
            plt.plot([frmaxs,frmaxs], np.log([min(s2n_spec), max(s2_spec)]), 'g')
            plt.text(10, max(np.log(s2_spec))-1, 'frmin = %.4f' % frmins)
            plt.text(10, max(np.log(s2_spec))-2, 'frmax = %.4f' % frmaxs)
            plt.xlim([0, int(max(s2_freq)+1)])
            plt.ylabel("ln(As) on" +chan[1]+' nm/s')
            sampling_rate = int(1/dt)
            plt.title('Station: %s, Sampleing rate: %s Hz' %(sta,sampling_rate))
            if franges2 >= franges1:
                plt.text(15, max(np.log(s1_spec))-1, 'prefered channel')
            plt.plot(freq_sx, np.log(spec_sx), 'k', label='signal in frange')
            try:
                plt.plot(freq_sx, syns, 'g--', linewidth=2, label='synthetic')
                plt.text(10, max(np.log(s1_spec))-3, 'slope = %.4f' % coeffs[0])
                plt.text(10, max(np.log(s1_spec))-4, 'residual = %.4f' % residus)
            except:
                pass
            plt.legend(loc='upper right')
            plt.savefig(snr_figure_dir+"/%s_%s_Ssnr.pdf" %(net,sta))
        
        if franges < snrcrtpara_s[1] or frmins > param['frmins']:
            goodS = False
            return goodP, goodS, spec, freq, n_spec, n_freq, spec_px, freq_px, spec_sx, freq_sx, frmin, frmax
        else:
            goodS = True
        
        if coeffs[0] < 0 and abs(residus) >= param['lincor'][1]:
            goodS = True
        else:
            goodS = False
            return goodP, goodS, spec, freq, n_spec, n_freq, spec_px, freq_px, spec_sx, freq_sx, frmin, frmax
    
    return goodP, goodS, spec, freq, n_spec, n_freq, spec_px, freq_px, spec_sx, freq_sx, frmin, frmax


def buildd(saving, stalst, ORIG, POS, icase, param, fc, allsite):
    """
    Build data matrix
        d = [ln(A1)-ln(C1)+ln(1+(f1i/fc)**2),                  ##
            ln(A2)-ln(C2)+ln(1+(f2i/fc)**2),                   ##
            ln(AM)-ln(CM)+ln(1+(fMi/fc)**2)]                   ##
    INPUT:  saving = saved spectrum for each station: saving[sta][1]['p']
            stalst = list of used stations
            fc     = corner frequency
            POS    = 'P' or 'S'
            icase  = 
                    1: high quality for finding best fc and alpha
                    2: low quality for t* inversion
                    3: low quality for t* inversion without bad fitting
            lnM    = when POS='S', log of seismic moment
    OUTPUT: data   = data matrix for t* inversion
    """
    
    if POS.upper() == 'P':
        ind = 0
    elif POS.upper() == 'S':
        ind = 1
    else:
        raise ValueError("P or S wave?")
    
    for ista in range(len(stalst)):
        sta = stalst[ista]
        freq_x = saving[sta][2][POS.lower()][0]
        spec_x = saving[sta][2][POS.lower()][1]
        if icase == 4:
            correc = saving[sta]['updated_corr'][ind]
            if param["add_site"] == 0:
                sitespec_x = 0
            elif sta in allsite[ind].keys():
                sitedata = allsite[ind][sta]
                sitefreq = sitedata[:, 0]
                sitespec = sitedata[:, 1]
                ind1 = np.nonzero(abs(freq_x[0]-sitefreq)<0.01)[0][0]
                ind2 = ind1 +len(freq_x)
                sitespec_x = sitespec[ind1:ind2]
            else:
                sitespec_x = 0
            
            stad = np.array([np.log(spec_x)+np.log(1+(freq_x/fc)**2)-sitespec_x-np.log(correc)]).transpose()
        else:  
            correc = saving[sta]['corr'][ind]
            stad = np.array([np.log(spec_x)+np.log(1+(freq_x/fc)**2)-np.log(correc)]).transpose()
              
        if ista == 0:
            data = stad
        else:
            data = np.vstack((data,stad))
    
    return data


def buildG(saving,stalst,alpha,POS,icase,param):

    if POS.upper() == 'P':
        ind = 0
    elif POS.upper() == 'S':
        ind = 1
    else:
        raise ValueError("P or S wave?")

    for ista in range(len(stalst)):
        sta = stalst[ista]
        freq_x = saving[sta][2][POS.lower()][0]
        exponent = -1*np.pi*freq_x*(freq_x**(-alpha))
        exponent = np.array([exponent]).transpose()
        Gblock = np.atleast_3d(exponent)

        if ista == 0:
            G = Gblock
        else:
            oldblock = np.hstack((G,np.zeros((G.shape[0], 1, 1))))
            newblock = np.hstack((np.zeros((Gblock.shape[0], G.shape[1], 1)),Gblock))
            G = np.vstack((oldblock,newblock))
    
    if icase == 2 or icase == 3:
        for i in range(G.shape[1]):
            row_num = np.count_nonzero(G[:,i,:])
            right_part = np.zeros((row_num,G.shape[1],1))
            right_part[:,i,:] = 1

            if i == 0:
                G_right_part = right_part
            else:
                G_right_part = np.vstack((G_right_part, right_part))
            
        G = np.hstack((G, G_right_part))

    if param['source_para'] == 1:    ## grid search for moment magnitude
        if POS.upper() == 'P':
            G = np.hstack((np.ones((G.shape[0], 1, 1)), G))
        else:
            G = np.hstack((np.ones((G.shape[0], 1, 1)), G))

    return G


def fitting(saving, sta, ORIG, POS, alpha, lnmomen, icase):
    """
    ## Calculate how well the synthetic spectrum fits the data
    """

    if POS.upper() == 'P':
        ind = 0
    elif POS.upper() == 'S':
        ind = 1
    else:
        raise ValueError("P or S wave?")
    
    if icase == 2 or icase == 3:
        corr = saving[sta]['corr'][ind]
        freq = saving[sta][2][POS.lower()][0]
        spec = saving[sta][2][POS.lower()][1]
        invtstar = saving[sta][icase]['tstar'][ind]
        inverror = saving[sta][icase]['error'][ind]
        newcorr = np.exp(np.log(corr)+inverror)
        synx = (newcorr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar)/(1+(freq/ORIG['fc'])**2))
        resid=(1-((np.linalg.norm(np.log(synx)-np.log(spec)))**2/(len(freq)-1)/np.var(np.log(spec)))) 
    
    elif icase == 4:
        corr = saving[sta]['updated_corr'][ind]
        freq = saving[sta][2][POS.lower()][0]
        spec = saving[sta][2][POS.lower()][1]
        invtstar = saving[sta][icase]['tstar'][ind]
        synx = (corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar)/(1+(freq/ORIG['fc'])**2))
        resid=(1-((np.linalg.norm(np.log(synx)-np.log(spec)))**2/(len(freq)-1)/np.var(np.log(spec)))) 
    
    return resid


def calresspec(saving, POS, lnmomen, fc, alpha):
    """
    ## Calculate residual spectrum for each station
    """

    if POS.upper() == 'P':
        ind = 0
    elif POS.upper() == 'S':
        ind = 1
    else:
        raise ValueError("P or S wave?")
    
    freq_x = saving[2][POS.lower()][0]
    spec_x = saving[2][POS.lower()][1]
    correc = saving['updated_corr'][ind]
    invtstar = saving[4]['tstar'][ind]
    righthand = lnmomen-np.pi*freq_x*(freq_x**(-alpha)*invtstar)
    resspec = np.array([np.log(spec_x)-np.log(correc)+np.log(1+(freq_x/fc)**2)-righthand])
    resratio = resspec/righthand*100
    resspec = np.vstack((freq_x, resspec))
    resspec = np.vstack((resspec, resratio))
    resspec = resspec.transpose()
    
    return resspec


def plotspec(saving, param, ARRIV, orid, POS, lnmomen, fc, alpha, icase, allsite):
    """
    Plot spectrum with good quality
    """

    net = ARRIV['net']
    sta = ARRIV['sta']
    if POS.upper() == 'P':
        ind = 0
        sitedata = allsite[ind]
    elif POS.upper() == 'S':
        ind = 1
        sitedata = allsite[ind]
    else:
        raise ValueError("P or S wave?")
    
    saving = saving[sta]
    spec = saving['spec'][ind]
    freq = saving['freq'][ind]
    n_spec = saving['nspec'][ind]
    n_freq = saving['nfreq'][ind]
    frmin = saving[2]['frmin'][ind]
    frmax = saving[2]['frmax'][ind]
    invtstar = saving[icase]['tstar'][ind]
    fitting = saving[icase]['fitting'][ind]
    updated_corr = saving['updated_corr'][ind]
    synspec = (updated_corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar)/(1+(freq/fc)**2))
    
    if POS.upper() == 'S':
        invtstarP = saving[icase]['tstar'][0]
        ttP = saving['Ptt']
        ttS = saving['Stt']
        QpQs = 2.25
        invtstar2 = invtstarP*QpQs*ttS/ttP
        synspec2 =(updated_corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar2)/(1+(freq/fc)**2))
        QpQs = 1.75
        invtstar2 = invtstarP*QpQs*ttS/ttP
        synspec3 = (updated_corr*np.exp(lnmomen)*np.exp(-np.pi*freq*(freq**(-alpha))*invtstar2)/(1+(freq/fc)**2))
    

    indx = np.all([(freq >= frmin), (freq < frmax)], axis=0)
    freqx = freq[indx]
    specx = spec[indx]
    synx = synspec[indx]
    resid = (1-((np.linalg.norm(np.log(synx)-np.log(specx)))**2/(len(freqx)-1)/np.var(np.log(specx)))) 

    # Plot spectrum for each station
    spec_figure_dir = tp.figdir3+"/%s" %(orid)
    if not os.path.exists(spec_figure_dir):
        os.makedirs(spec_figure_dir)
    fig = plt.figure(7, figsize=(8, 6))
    plt.clf()
    ax = fig.add_subplot(111)
    if param['plotspecloglog']:
        plt.loglog(freq, spec, 'b', label='signal')
        plt.loglog(n_freq, n_spec, 'r', label='noise')
        plt.loglog([frmin, frmin], [min(n_spec), max(spec)], 'g', label='frmin')
        plt.loglog([frmax, frmax], [min(n_spec), max(spec)], 'g', label='frmax')
        plt.loglog(freq, synspec, 'g--', linewidth=2, label='synthetic')
        plt.legend(loc='lower left')
        text =  "t* = %.2f\n" %(invtstar) +\
                "frange = [%.2f, %.2f]\n" %(frmin, frmax)+ \
                "sd = %.4f\n" %(saving[icase]['err'][ind]) + \
                "fitting = %.4f\n" %(fitting)

        plt.text(0.6, 0.75, text, transform=ax.transAxes)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('log(A%s), nm/s' %(POS.lower()))
        plt.title('%s Station: %s' % (orid, sta))
        plt.savefig(spec_figure_dir+'/%s_%s_%sspectrum_loglog.pdf' % (net, sta, POS.upper()))
    else:
        if param["add_site"] == 0:
            sitespec_x = 0
        elif sta in sitedata.keys():
            sitefreq = sitedata[sta][:, 0]
            sitespec = sitedata[sta][:, 1]
            ind1 = np.nonzero(abs(freqx[0]-sitefreq)<0.01)[0][0]
            ind2 = ind1+len(freqx)
            sitespec_x = sitespec[ind1:ind2]
            plt.plot(freqx, np.log(specx)-sitespec_x,'k--',label='correctsite_signal')
        else:
            sitespec_x = 0
        
        plt.plot(freq, np.log(spec), 'b', label='signal')
        plt.plot(n_freq, np.log(n_spec), 'r', label='noise')
        plt.plot([frmin, frmin], np.log([min(n_spec), max(spec)]), 'g', label='frmin')
        plt.plot([frmax, frmax], np.log([min(n_spec), max(spec)]), 'g', label='frmax')
        plt.plot(freq, np.log(synspec), 'g--', linewidth=2, label='synthetic')
        if POS.upper() == "S":
            invtstarP = saving[icase]['tstar'][0]
            ttP = saving['Ptt']
            ttS = saving['Stt']
            QpQs = 2.25
            invtstar2 = invtstarP*QpQs*ttS/ttP
            synspec2 =(updated_corr*np.exp(lnmomen)*np.exp(-np.pi*freqx*(freqx**(-alpha))*invtstar2)/(1+(freqx/fc)**2))
            QpQs = 1.75
            invtstar2 = invtstarP*QpQs*ttS/ttP
            synspec3 = (updated_corr*np.exp(lnmomen)*np.exp(-np.pi*freqx*(freqx**(-alpha))*invtstar2)/(1+(freqx/fc)**2))
            plt.plot(freqx, np.log(synspec2), 'm--', linewidth=2, label='QpQs=2.25')
            plt.plot(freqx, np.log(synspec3), 'c--', linewidth=2, label='QpQs=1.75')
        
        plt.legend(loc='lower left')
        text =  "t* = %.2f\n" %(invtstar) +\
                "frange = [%.2f, %.2f]\n" %(frmin, frmax)+ \
                "sd = %.4f\n" %(saving[icase]['err'][ind]) + \
                "fitting = %.4f\n" %(fitting) +\
                "redidusl = %.4f\n" %(resid)

        plt.text(0.6, 0.75, text, transform=ax.transAxes)
        plt.text(0.6, 0.75, text, transform=ax.transAxes)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('log(A%s), nm/s' %(POS.lower()))
        plt.title('%s  Station: %s.%s' % (orid, net, sta))
        plt.savefig(spec_figure_dir+'/%s_%s_%sspectrum_loglinear.pdf' % (net, sta, POS.upper()))
        
        # Output for plotting
        if POS.upper() == "P":
            spec_file = sta+"_signal_spec.dat"
            n_spec_file = sta+"_noise_spec.dat"
            syn_spec_file = sta+"_syn_spec.dat"
            fspec = open(spec_file, "w")
            fnspec = open(n_spec_file, "w")
            fsynspec = open(syn_spec_file, "w")
            for i in range(len(freq)):
                string1 = str(freq[i])+" "+str(np.log(spec)[i])
                fspec.write(string1+"\n")
                string2 = str(freq[i])+" "+str(np.log(n_spec)[i])
                fnspec.write(string2+"\n")
                string3 = str(freq[i])+" "+str(np.log(synspec)[i])
                fsynspec.write(string3+"\n")
            
            fspec.close()
            fnspec.close()
            fsynspec.close()