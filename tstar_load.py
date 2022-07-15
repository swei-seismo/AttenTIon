#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import subprocess as sp
import numpy as np
import tstar_parameters as tp

def channels(sta):
    """
    Determine the channel names for each station, need to modify for your dataset
    OUTPUT: (chan) = Chanels of P and S waves
    """
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
    return chan

def loadORIG(orid, source_para,ORIG,cmd):
    """
    :param source_para:
        The option for determining the source spectra-
            1: grid search for fc and Mw, need to load mb from './data/eventid_sub.ls'
            2: need to load fc and Mw from './data/TongaEQ_source_spectra_v4.dat'
            3: input source spectra directly. No idea of how to deal with it now
    OUTPUT: (ORIG) = Origin Disctionary with Keys:
                [orid] = origin id (orid)
                [time] = origin time (unix seconds)
                [lat]  = hypocenter latitude
                [lon]  = hypocenter longitude
                [dep]  = hypocenter depth
                [mb]   = local magnitude
                [mw]   = seismic magnitude
                [fc]   = corner frequency
    """
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)      
    output = p.communicate().__str__()
    (evla, evlo, evdp) = (float(output.split()[1]), float(output.split()[2]), \
        output.split()[3])
    # evdp = output.split()[9]        ##UNITS: km
    evdp = float(evdp[:evdp.find('\\n')]) 
    # evdp = float(evdp)/1000     ##UNITS: km

    (ORIG['orid'],ORIG['lat'],ORIG['lon'],ORIG['dep'],ORIG['mb'],ORIG['mw'],ORIG['mo'],ORIG['fc'])\
        =(orid, evla, evlo, evdp, -1, -1, -1, -1 )

    ######need to modify the input files and formats for your dataset######
    if source_para == 1:
        sub = orid.split('_')[0]
        evlst = os.path.abspath(os.path.dirname(os.path.dirname(os.getcwd())))+"/eventid_sub%s.lst"%sub
        for line in open(evlst).readlines()[1:]:
            if orid.split('_')[1] == line.split()[0]:
                ORIG['mb'] = float(line.split()[6])
                break
        open(evlst).seek(0)
    elif source_para == 2:
        SSlst = os.path.abspath(os.path.dirname(os.path.dirname(os.getcwd())))+"/TongaEQ_source_spectra_v4.dat"
        for line in open(SSlst).readlines():
            if orid == line.split()[0]:
                ORIG['mw'] = float(line.split()[4])
                ORIG['mo'] = pow(10, 1.5*ORIG['mw']+9.095)
                ORIG['fc'] = float(line.split()[5])
                break
        open(SSlst).seek(0)
    # elif source_para == 3:
    #     help!!
    ########################################################################
    return ORIG

def loaddata(param, orid, fldir):
    """
    Load data from 
    1) sac files from './data/processedSeismograms'
    2) event body wave magnitude 
    and *optionally 
        1) load fc and Mw from input file './data/' 
        2) load P and S arrival times from './data/arriv.dat'

    INPUT: 
       (orid) = origin id (orid)
       (fldir) = sac file directory

    OUTPUT: stalst = station list for this event
            (ORIG) = Origin Disctionary 
            (ARRIV) = Arrival Dictionary with Keys:
                [sta]   = station name
                [ptime] = P arrival time (unix seconds)--to origin
                [parid] = P arrival id (arid)
                [pchan] = channel that P was picked on
                [stime] = S arrival time (unix seconds)--to origin
                [sarid] = S arrival id (arid)
                [schan] = channel that S was picked on
                [delta] = distance in degrees from sta to origin
                [baz]   = station to event azimuth
            (Ptt) = P travel time
            (Stt) = S travel time
            (SDATA) = If S arrival exists
    
    """
    fls = os.listdir(fldir)
    ##ORIGIN INFORMATION
    ORIG = {}
    ARRIV = {}

    stalst = [ista.split('_')[1].split('.')[0] for ista in fls]
    stalst = list(set(stalst))
    if len(stalst) == 0:
        return 0, ORIG, ARRIV
    
    os.chdir(fldir)
    for i, sta in enumerate(stalst):
        chan = channels(sta)
        ##load event and station information
        if i == 0: 
            ##load event information in the first station
            cmd = 'saclst evla evlo evdp f *%s*%s*sac'%(sta,chan[0])
            ORIG = loadORIG(orid,param['source_para'],ORIG,cmd)
        cmd = 'saclst T0 T1 dist baz stla stlo f *%s*%s*sac'%(sta,chan[0])
        p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)      
        output = p.communicate().__str__()
        (T0, T1, dist, baz, stla, stlo) = (output.split()[1], output.split()[2], output.split()[3], \
            output.split()[4], output.split()[5], output.split()[6])
        if T0 == '-12345':
            T0 = -1
        if T1 == '-12345':
            T1 = -1

        ######need to modify the input file and data formats for your dataset######
        if param['input_arriv'] == 1:
            sub = orid.split('_')[0]
            arrivlst = os.path.abspath(os.path.dirname(os.path.dirname(os.getcwd())))+"/PSass_%s.lst"%sub
            for line in open(arrivlst).readlines()[1:]:
                if orid.split('_')[1] == line.split()[0] and sta == line.split()[5]:
                    (ortime, ptime, stime) = \
                        (float(line.split()[4]), float(line.split()[9]), float(line.split()[11]))
                    T0 = (ptime == -1 and -1 or ptime-ortime)
                    T1 = (stime == -1 and -1 or stime-ortime)
                    break
            open(arrivlst).seek(0)
        ###########################################################################
        if T1 == -1 :
            SDATA = False
        else:
            SDATA = True
        ARRIV.setdefault(sta,{})['T0'] = float(T0)
        ARRIV.setdefault(sta,{})['T1'] = float(T1)
        ARRIV.setdefault(sta,{})['delta'] = 360*float(dist)/(np.pi*2*param["Re"])
        ARRIV.setdefault(sta,{})['baz'] = float(baz)
        ARRIV.setdefault(sta,{})['pchan'] = chan[2]
        ARRIV.setdefault(sta,{})['SDATA'] = SDATA

        # print(ARRIV[sta]['SDATA'],ARRIV[sta]['T1'])
    numSTA = len(stalst)
    numP = len([sta for sta in stalst if ARRIV[sta]['T0']!=-1])
    numS = len([sta for sta in stalst if ARRIV[sta]['T1']!=-1])
    if numSTA > numS:
        tp.logfl.write('%d Records: %d for P and %d for S\n' % (numSTA,numP,numS))
    elif numSTA < numS:
        tp.logfl.write('%d Records: %d for P and %d for S\n' % (numSTA,numP,numS))
        tp.logfl.write('You must have numP >= numS\n')
        stalst = 0

    

    return stalst, ORIG, ARRIV

def loadGS(orid, gsdir):
    """
    LOAD FILES WITH GS (GEOMETRICAL SPREADING) AND FREE SURFACE EFFECT,
       ETC FOR ZERO FREQ. AMPLITUDE
    """
    pgsfile='%s/pgsfile%s.txt' % (gsdir,orid)
    if not os.path.isfile(pgsfile):
        print('File %s does not exist, using mkgsfl.sh to generate' % pgsfile)
        exit()
    PGS={'gval':[line.split()[0] for line in open(pgsfile)],'stalist':[line.split()[1] for line in open(pgsfile)]}
    
    sgsfile='%s/sgsfile%s.txt' % (gsdir,orid)
    if not os.path.isfile(sgsfile):
        print('File %s does not exist, using mkgsfl.sh to generate' % sgsfile)
        exit()
    SGS={'gval':[line.split()[0] for line in open(sgsfile)],'stalist':[line.split()[1] for line in open(sgsfile)]}

    return PGS, SGS