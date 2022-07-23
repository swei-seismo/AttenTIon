#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import glob
import tstar_parameters as tp


def loadchannels(orid, net, sta):
    """
    Determine the channel names of each station
    INPUT:  (orid) ==> Event ID
            (net)  ==> Network name
            (sta)  ==> Station name
    OUTPUT: (chan) ==> A list contain channel names
    """
    sacfiles = tp.sacdir+"/%s/%s.%s.*.sac" %(orid,net,sta)
    command = "saclst KCMPNM f %s | awk -F ' ' '{print $2}'" %(sacfiles)
    output = os.popen(command).read().strip().split('\n')
    chan = sorted(output) # confirm chan[0]-E, chan[1]-N, chan[2]-Z

    return chan


def loadORIG(orid, param):
    """
    Load event info of each event
    :param source_para:
        The option for determining the source spectra-
            1: grid search for fc and Mw, need to load mb from './data/eventid_sub.ls'
               (need to specify the path to this file in this function)
            2: need to load fc and Mw from './data/TongaEQ_source_spectra_v4.dat'
               (need to specify the path to this file in this function) 
            3: input source spectra directly. No idea of how to deal with it now
    
    OUTPUT: (ORIG) = Origin Disctionary with Keys:
                [orid] = origin id (orid)
                [time] = origin time (unix seconds) # Obspy.UTCDateTime?
                [lat]  = hypocenter latitude
                [lon]  = hypocenter longitude
                [dep]  = hypocenter depth
                [mb]   = local magnitude
                [mw]   = seismic magnitude
                [fc]   = corner frequency
    """
    ## A loop may be used.
    ORIG = {}
    catalog_file = "./data/eventid.lst"
    command = "cat %s | awk '/%s/ {print}'" %(catalog_file, orid)
    event_info = os.popen(command).read().strip()
    # Need to modify the colunm number
    evlo = float(event_info.split()[2])
    evla = float(event_info.split()[1])
    evdp = float(event_info.split()[3]) ## Unit: km

    ORIG['orid'] = orid
    ORIG['lat'] = evla
    ORIG['lon'] = evlo
    ORIG['dep'] = evdp
    ORIG['mb'] = -1
    ORIG['mw'] = -1
    ORIG['mo'] = -1
    ORIG['fc'] = -1

    ###### Need to modify the input files and formats for your dataset ######
    if param['source_para'] == 1:
        ORIG['mb'] = float(event_info.split()[6])
    
    elif param['source_para'] == 2:
        SSlst = "./data/TongaEQ_source_spectra_v4.dat"
        cmd = "cat %s | awk '/%s/ {print}'" %(SSlst, orid)
        output = os.popen(cmd).read().strip()
        ORIG['mw'] = float(output.split()[4])
        ORIG['mo'] = pow(10, 1.5*ORIG['mw']+9.095)
        ORIG['fc'] = float(output.split()[5])
    
    # elif source_para == 3:
    #     help!!
    #########################################################################
    
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
       (param) = parameters in tstar_parameters.py 
       (orid) = origin id (orid)
       (fldir) = sac file directory ('./data/processedSeismograms/orid')

    OUTPUT: (stalst) = station list for this event 
            (ARRIV)  = Arrival Dictionary with Keys:
                [net]   = network name
                [sta]   = station name
                [T0]    = P arrival time
                [T1]    = S arrival time
                [dist]  = distace in km from station to event
                [baz]   = station to event azimuth ? (event to station or station to event?)
                [gcarc] = distance in degree from station to event
                [dt]    = sampling rate
                [chan]  = channel names for each station
                [pchan] = channel that P arrival was picked on
                [SDATA] = If S arrival exists
            
            An example for ARRIV:
                ARRIV = {'EP15':{'T0':10,
                                 'T1':15},
                         'KD02':{'T0':20,
                                 'T1':25}
                        }
                #[parid] = P arrival id (arid)
                #[sarid] = S arrival id (arid)
    """
    os.chdir(fldir)
    sacfls = glob.glob("*Z.sac")
    stalst = [sacfl.split(".")[1] for sacfl in sacfls]
    ARRIV = {}
    
    if len(stalst) == 0:
        return stalst, ARRIV
    
    for sacfl in sacfls:
        net = sacfl.split(".")[0]
        sta = sacfl.split(".")[1]
        chan = loadchannels(orid, net, sta)
        cmd = "saclst t0 t1 dist baz gcarc delta f %s" %(sacfl)
        output = os.popen(cmd).read().strip()
        T0 = output.split()[1]      
        T1 = output.split()[2]    
        dist = output.split()[3]
        baz = output.split()[4]
        gcarc = output.split()[5]
        dt = output.split()[6]
        if T0 == '-12345':
            T0 = -1
        if T1 == '-12345':
            T1 = -1
        
        ######need to modify the input file and data formats for your dataset######
        ## We can also use 'awk' command to find the corresponding values.
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
        
        if T1 == -1:
            SDATA = False
        else:
            SDATA = True
        ARRIV.setdefault(sta,{})['net'] = net
        ARRIV.setdefault(sta,{})['sta'] = sta
        ARRIV.setdefault(sta,{})['T0'] = float(T0)
        ARRIV.setdefault(sta,{})['T1'] = float(T1)
        ARRIV.setdefault(sta,{})['dist'] = float(dist)
        ARRIV.setdefault(sta,{})['baz'] = float(baz)
        ARRIV.setdefault(sta,{})['gcarc'] = float(gcarc)
        ARRIV.setdefault(sta,{})['dt'] = float(dt)
        ARRIV.setdefault(sta,{})['chan'] = chan
        ARRIV.setdefault(sta,{})['pchan'] = chan[2]
        ARRIV.setdefault(sta,{})['SDATA'] = SDATA

    numSTA = len(stalst)
    numP = len([sta for sta in stalst if ARRIV[sta]['T0'] != -1])
    numS = len([sta for sta in stalst if ARRIV[sta]['T1'] != -1])
    if numSTA >= numS:
        tp.logfl.write('%d Records: %d for P and %d for S\n' % (numSTA,numP,numS))
    elif numSTA < numS:
        tp.logfl.write('%d Records: %d for P and %d for S\n' % (numSTA,numP,numS))
        tp.logfl.write('You must have numP >= numS\n')
        stalst = []

    return stalst, ARRIV


def loadGS(orid, gsdir):
    """
    Load files with GS (Geometrical Spreading) and free surface effect,
        etc for zero freq. Amplitude ??
    """
    pgsfile = '%s/pgsfile%s.txt' % (gsdir,orid)
    if not os.path.isfile(pgsfile):
        print('File %s does not exist, using mkgsfl.py to generate' % pgsfile)
        exit()
    PGS = {'gval':[line.split()[0] for line in open(pgsfile)],
           'stalist':[line.split()[1] for line in open(pgsfile)]}
    
    sgsfile = '%s/sgsfile%s.txt' % (gsdir,orid)
    if not os.path.isfile(sgsfile):
        print('File %s does not exist, using mkgsfl.py to generate' % sgsfile)
        exit()
    SGS = {'gval':[line.split()[0] for line in open(sgsfile)],
           'stalist':[line.split()[1] for line in open(sgsfile)]}

    return PGS, SGS