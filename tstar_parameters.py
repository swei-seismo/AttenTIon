#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Set parameters used in the inversion 

import os

def set_parameters():
    """
    Parameters used in frequency and time domains, modify them to fit your dataset.

    :param snrcrtp1: Only matters when source_para = 1.

    :param dstress:
        Stress drop in [MPa], for calculating the reasonable range of theoretical corner frequencies.

    :param minf and maxf:
        Find the longest segment that meet the SNR thresold.
    
    :param mtspec_para:
        The option for determining which function will be used in multitaper
        1: mtspec
        2: sine_psd
    
    :param doplotseis: bool
        Determine whether to plot all the three-component seismograms in ./workdir/specfig/plotseis.
    :param doplotsnr: bool
        Determine whether to plot all the SNR of spectra in ./workdir/specfig/plotsnr.
    :param doplotspec: bool
        Determine whether to plot all the spectra in log(Am)-linear(Fre) scale in ./workdir/specfig/plotspec.
    :param plotspecloglog: bool
        Determine whether to plot all the spectra in log(Am)-log(Fre) scale in ./workdir/specfig/plotspec.
    :param doplotfcts: bool
        Determine whether to plot t* perturbation versus corner frequency in ./workdir/specfig/plottstar-fc.
    :param doplotfcall: bool
        Determine whether to plot L2 Norm and Seismic Moment versus corner frequency in ./workdir/specfig/plotfcall.

    :param source_para:
        The option for determining the source spectra-
            1: grid search for fc and Mw to calculate the source spectra.
            2: input fc and Mw to caulcte the source spectra.
            3: input source spectra directly.

    :param input_arriv:
        The option to read P and S arrival times-
            0: read arrival time from sac files
            1: read arrival time from input file './data/arriv.dat'
   """
    param = {
        "snrcrtp1": [3.0, 3.0, 3.0],  # Minimum [SNR, frequency band width] for finding best fc and alpha
                                      # See the hints for the second and third parameter in tstarsub.longseg
        "snrcrtp2": [3.0, 3.0, 3.0],  # Minimum [SNR, frequency band width] for t* inversion
        "snrcrts1": [2.0, 2.0, 2.0],
        "snrcrts2": [2.0, 2.0, 2.0],
        "lincor": [0.7, 0.7, 0.7],    # Minimum linear correlation coefficients
        "misfitP": 0.75,              # Minimum misfit for P wave measured by correlative function
        "misfitP2": 0.2,              # Maximum L2 Norm for P wave
        "misfitS": 0.75,              # Minimum misfit for S wave measured by correlative function
        "misfitS2": 0.2,              # Maximum L2 Norm for S wave

        "dstress": [0.1, 100.0],

        # Time domain parameters
        "WLP": 3.0,         # Window length of P wave signal and corresponding noise
        "WLS": 6.0,         # Window length of S wave signal and corresponding noise
        "prewin": [0.5,1],  # Window start before arrival
        "gaps": 0.5,        # Seconds between S noise and signal
        "beta": 4000,       # Shear velocity in m/s for approximating corner frequency

        # Frequency domain parameters
        "frminp": 8,
        "frmins": 8,
        "40": [0.05, 16.05],
        "50": [0.05, 20.05],
        "100": [0.05, 40.05],

        # Options for multitaper method
        "mtspec_para":2, 

        # Options for plotting figures

        "doplotseis": False,
        "doplotsnr": False,
        "doplotspec": True,
        "plotspecloglog": False,
        "doplotfcts": False,
        "doplotfcall": False,

        # "doplotseis": False,
        # "doplotsnr": True,
        # "doplotspec": True,
        # "plotspecloglog": False,
        # "doplotfcts": False,
        # "doplotfcall": False,

        # Options for determining source parameters
        "source_para": 1,
        "input_arriv": 0, 
        "alpha": 0.27,
        "fcps": 1.0,       ## fc(P)/fc(S)

        # Options for correcting site effects (1-correct, 0-no correction)
        "add_site": 0

    }

    return param


def working_dir():
    """ 
    Define all the working dictionaries with absolute paths and return the event list
    """

    global workdir, sacdir, gsdir, figdir, resultdir, sitedir, figdir1, \
           figdir2, figdir3, figdir4, figdir5, logfl, fclist, catalog, spectra, fittingfl
    
    maindir = "/mnt/research/seismo_wei/ZZhang/tstar_inversion"
    sacdir = '/mnt/research/seismo_wei/ZZhang/tstar_inversion/tstar_data/correctedSeismograms'
    namedir = 'Eventname'
    gsdir = maindir + '/tstar_data/GS'
    catalog = maindir + '/tstar_data/AACSE_catalog.dat'
    spectra = maindir + '/tstar_data/Alaska_events_source_spectra.dat'
    
    ## Output
    currentdir = os.getcwd()
    sitedir = currentdir+"/data/sitespec/"
    workdir = currentdir+'/workdir' ## Output dictionary
    if not os.path.exists(workdir):
        os.makedirs(workdir, exist_ok=True)
    
    resultdir = workdir + '/result'
    sitedir = resultdir + '/site' 
    figdir  = workdir + '/specfig'
    figdir1 = figdir + '/plotseis'
    figdir2 = figdir + '/plotsnr'
    figdir3 = figdir + '/plotspec' 
    figdir4 = figdir + '/plottstar-fc'
    figdir5 = figdir + '/plotfall'

    dir_lst = [figdir, resultdir, sitedir, figdir1, figdir2, figdir3, figdir4, figdir5]
    for idir in dir_lst:
        if not os.path.exists(idir):
            os.makedirs(idir, exist_ok=True)

    logfile = workdir + '/event.log' 
    logfl = open(logfile,'a')
    fclist = open(workdir + '/bestfc.lst','a')
    fittingfile = workdir + '/fitting.lst'
    fittingfl = open(fittingfile, 'a') 
    
    oridlst = []
    forid = open(namedir, 'r')
    for name in forid.readlines():
        name = name.strip()
        if name is not None:
            oridlst.append(name)
    forid.close()

    return oridlst
