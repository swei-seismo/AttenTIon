#!/usr/bin/env python3
# -*- coding:utf-8 -*-
## SET PARAMETERS USED IN THE INVERSION 

import os

def set_parameters():
    """
    Parameters used in frequency and time domains, modify them to fit your dataset.

    :param snrcrtp1: Only matters when source_para = 1.
        Minimum [SNR, frequency band width???] for finding best fc and alpha.

    :param dstress:
        Stress drop in [MPa], for calculating the reasonal range of theoretical corner frequencies.

    :param doplotseis: bool
        Determine whether to plot all the three-component seismograms in ./wordir/specfig/plotseis.
    :param doplotsnr: bool

    :param source_para:
        The option for determining the source spectra-
            1: grid search for fc and Mw to calculate the source spectra
            2: input fc and Mw to caulcte the source spectra
            3: input source spectra directly.
    :param input_arriv:
        The option to read P and S arrival times-
            0: read arrival time from sac files
            1: read arrival time from input file './data/arriv.dat'
   """
    param = {
        ## FREQUENCY DOMAIN PARAMETERS AND ALPHA
        # "loco": 0.05,
        # "hico": 9.9,
        # "npol": 4,
        "frminp": 4,

        "snrcrtp1": [3.0,4.0,3.0],   ## MINIMUM [SNR,FREQUENCY BAND WIDTH] FOR FINDING BEST fc AND alpha
        ## see the hints for the second and third parameter in longseg
        "snrcrtp2": [2.0,3.0,2.0],    ## MINIMUM [SNR,FREQUENCY BAND WIDTH] FOR t* INVERSION
        "snrcrts1": [2.5,2.0,1.2],    ## MINIMUM [SNR,FREQUENCY BAND WIDTH] FOR FINDING BEST fc AND alpha
        "snrcrts2": [1.8,1.0,1.2],    ## MINIMUM [SNR,FREQUENCY BAND WIDTH] FOR t* INVERSION
        "lincor": [0.7,0.7,0.7],    ## MINIMUM LINEAR CORRELATION COEFFICIENTS
        "misfitP": 0.85,    ## MINIMUM MISFIT FOR P WAVE MEASURED BY CORRELATIVE FUNCTION
        "misfitP2": 0.2,    ## MAXIMUM L2 NORM FOR P WAVE
        "misfitS": 0.85,    ## MINIMUM MISFIT FOR S WAVE MEASURED BY CORRELATIVE FUNCTION
        "misfitS2": 0.2,    ## MAXIMUM L2 NORM FOR S WAVE

        "dstress": [0.5,20.0],   
        "Re": 6371,

        ##TIME DOMAIN PARAMETERS
        "WLP": 5.0,
        "WLS": 8.0,
        "pretime": 50,
        "posttime": 150,
        "prewin": [0.5,1],  ## WINDOW START BEFORE ARRIVAL
        "gaps": 0.5,        ## seconds between S noise and signal
        "dt": 0.025,
        "beta": 4000,       ## SHEAR VELOCITY IN m/s FOR APPROXIMATING CORNER FREQUENCY

        ##OPTIONS FOR PLOTTING FIGURES
        "doplotseis": False,
        "doplotsnr": False,
        "doplotspec": False,
        "plotspecloglog": True,
        "doplotfcts": True,
        "doplotfcall": True,

        ##OPTIONS FOR DETERMINING SOURCE PARAMTERS
        "source_para": 2,
        "input_arriv": 1, 
        "alpha": 0.27,
        "fcps": 1.0,       ## fc(P)/fc(S)

    }
    return param

def working_dir():
    """ All the working directories with absolute paths
    """
    global workdir, sacdir, gsdir, figdir, resultdir, figdir1, figdir2, figdir3, figdir4, figdir5, logfl, fclist
    maindir = os.path.abspath(os.getcwd())
    # workdir = maindir + '/workdir'    ## OUTPUT DIRECTOR
    
    # need to modify! (please use: workdir = maindir + '/workdir')
    workdir = '/mnt/gs21/scratch/yurong/tstar/workdir' # OUTPUT to scratch dorectory on hpcc when running with a large dataset
    # need to modify! (please use: workdir = maindir + '/workdir')
    if not os.path.isdir(workdir):
        os.mkdir(workdir)   

    namedir = maindir + '/data/Eventname'
    # sacdir = maindir + '/data/processedSeismograms'
    # need to modify! (please use: sacdir = maindir + '/data/processedSeismograms')
    sacdir = '/mnt/home/yurong/Research/tstar/4YZhang/sacfl'
    # need to modify! (please use: sacdir = maindir + '/data/processedSeismograms')

    gsdir = maindir + '/data/GS'
    
    resultdir = workdir + '/result' 
    figdir  = workdir + '/specfig'
    figdir1 = figdir + '/plotseis'
    figdir2 = figdir + '/plotsnr'
    figdir3 = figdir + '/plotspec' 
    figdir4 = figdir + '/plottstar-fc'
    figdir5 = figdir + '/plotfall'

    logfile = workdir + '/event.log' 
    logfl=open(logfile,'a')
    fclist=open(workdir + '/bestfc.lst','a')
    

    dir_lst = [figdir, resultdir, figdir1, figdir2, figdir3, figdir4, figdir5]
    for idir in dir_lst:
        if not os.path.isdir(idir):
            os.mkdir(idir)
    
    oridlst=[]
    for name in open(namedir).readlines():
        name = name.strip('\n')
        if name is not None:
            oridlst.append(name)

    return oridlst

