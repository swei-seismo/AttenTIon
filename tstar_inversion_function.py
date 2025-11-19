#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Main function for t* inversion 

import os
import tstarsub
import numpy as np
from scipy import optimize
import tstar_parameters as tp


def Loop_Record(orid, stalst, param, ARRIV, PGS, SGS, ORIG):
    """
    Loop over each record, process the seismic information
    """

    saving = {}       # Saved spectra: saving[sta][icase][POS]
    staP2lst = []     # Stations used in t*(P) inversion
    staP3lst = []     # Stations used in t*(P) inversion without bad fitting
    staS2lst = []     # Stations used in t*(S) inversion
    staS3lst = []     # Stations used in t*(S) inversion without bad fitting
    for i, sta in enumerate(stalst):

        # print('Working on station %s (%d of %d)' % (sta, i+1, len(stalst)))
        chan = ARRIV[sta]['chan']

        # Return raw data (P arrival is required, but S arrival is not necessary)
        (dd, tt, flag) = tstarsub.readseismo(orid, ARRIV[sta])
        if not flag:
            print('ERROR: Unable to read %s.%s' % (orid,sta))
            continue

        # Determine which horizontal channel of S arrival has high signal-to-noise ratio
        if ARRIV[sta]['SDATA']:
            stiim = ARRIV[sta]['T1']-ARRIV[sta]['T0']   # S arrival in dd series
            indn = np.nonzero((tt>(stiim-1-param['WLS']))&(tt<stiim-1))
            inds = np.nonzero((tt>stiim)&(tt<(stiim+param['WLS'])))
            if np.absolute(dd[1][inds]).max() == 0:
                ARRIV.setdefault(sta,{})['schan'] = chan[0]
            elif np.absolute(dd[0][inds]).max() == 0:
                ARRIV.setdefault(sta,{})['schan'] = chan[1]
            else:
                snrew = np.absolute(dd[0][inds]).max()/np.absolute(dd[0][indn]).max()
                snrns = np.absolute(dd[1][inds]).max()/np.absolute(dd[1][indn]).max()
                if snrew > snrns:
                    ARRIV.setdefault(sta,{})['schan'] = chan[0]
                else:
                    ARRIV.setdefault(sta,{})['schan'] = chan[1]

        (p_dd, pn_dd, s_dd1, sn_dd1, pc_dd1, s_dd2, sn_dd2, pc_dd2) = tstarsub.fixwin(dd, tt, param, ARRIV[sta], orid)
        PWINDATA = np.vstack((p_dd, pn_dd))
        SWINDATA1 = np.vstack((s_dd1, sn_dd1))
        SWINDATA1 = np.vstack((SWINDATA1, pc_dd1))
        SWINDATA2 = np.vstack((s_dd2, sn_dd2))
        SWINDATA2 = np.vstack((SWINDATA2, pc_dd2))

        # Calculate spectra and auto selects frequency band above set SNR
        ##======== 2 means lower quality data for t* inversion ========##
        (goodP2, goodS2, spec, freq, n_spec, n_freq, spec_px, freq_px, spec_sx, freq_sx, frmin, frmax) \
            = tstarsub.dospec(PWINDATA, SWINDATA1, SWINDATA2, orid, param, 2, ARRIV[sta], ORIG)
        if not goodP2:
            # print('No good P wave signal. Skip to next record.')
            continue

        # Save spectrum and other information for each station
        saving[sta] = {}
        saving[sta]['spec'] = spec
        saving[sta]['nspec'] = n_spec
        saving[sta]['freq'] = freq
        saving[sta]['nfreq'] = n_freq
        saving[sta][1], saving[sta][2], saving[sta][3], saving[sta][4] = {}, {}, {}, {}

        saving[sta]['Ptt'] = ARRIV[sta]['T0']-ARRIV[sta]['otime']
        saving[sta]['Stt'] = ARRIV[sta]['T1']-ARRIV[sta]['otime']
        saving[sta][2]['good'] = [goodP2, goodS2]
        saving[sta][2]['frmin'], saving[sta][2]['frmax'] = frmin, frmax
        saving[sta][2]['p'] = [freq_px, spec_px]
        saving[sta][2]['s'] = [freq_sx, spec_sx]

        # Corrections of GS
        correcP = float(PGS['gval'][PGS['stalist'].index(sta)])
        if correcP == 0:
            print('Bad correction of P wave for station %s' % sta)
            continue
        saving[sta]['corr'] = [correcP]
        staP2lst.append(sta)

        if ARRIV[sta]['SDATA'] and goodS2:
            correcS = float(SGS['gval'][SGS['stalist'].index(sta)])
            if correcS == 0:
                print('Bad correction of S wave for station %s' % sta)
                continue
            saving[sta]['corr'].append(correcS)
            staS2lst.append(sta)
        ##======== 2 means lower quality data for t* inversion ========##

    return staP2lst, staP3lst, staS2lst, staS3lst, saving, ARRIV


def fixed_bestfc(orid, ORIG):
    """
    If corner frequency is given, like from decomposition method, then use the fixed corner frequency.
    """

    cmd = "awk '/%s/ {print $6}' %s" %(orid, tp.spectra)
    output = os.popen(cmd).read().strip()
    ORIG['fc'] = float(output)
    # ORIG['fc'] = 100

    return ORIG


def get_good_stalst(orid, saving, stalst, ORIG, POS, icase, param, allsite):
    """
    #################### INVERSION FOR t* ############################
    ##    EQUATION 3 IN Stachnik, Abers, Christensen, 2004, JGR     ##
    ##      d = Gm      (Nx1) = (Nx(M+1+)) ((M+1)x1)                ##
    ##      FOR A GIVEN fc AND alpha (GRID SEARCH):                 ##
    ##      d = [ln(A1)-ln(C1)+ln(1+(f1i/fc)**2),                   ##
    ##           ln(A2)-ln(C2)+ln(1+(f2i/fc)**2),                   ##
    ##           ln(AM)-ln(CM)+ln(1+(fMi/fc)**2)]                   ##
    ##      G = [[1, -pi*f1i*f1i**(-alpha), 0, ..., 0],             ##
    ##           [1, 0, -pi*f2i*f2i**(-alpha), ..., 0],             ##
    ##           [1, 0, 0, ..., -pi*fMi*fMi**(-alpha)]]             ##
    ##      m = [[ln(Moment)],[tstar01],[tstar02],...,[tstar0M]     ##
    ##################################################################
    """

    ## icase = 2
    if POS == "P":
        bestfc = ORIG['fc']
    elif POS == "S":
        bestfc = ORIG['fc']/param['fcps']
    
    data = tstarsub.buildd(saving, stalst, ORIG, POS, icase, param, bestfc, allsite)
    G = tstarsub.buildG(saving, stalst, param['alpha'], POS, icase, param)
    
    l_bound = [-np.inf]*(2*len(stalst)+1)
    l_bound[0] = 21
    u_bound = [np.inf]*(2*len(stalst)+1)
    u_bound[0] = 40
    result = optimize.lsq_linear(G[:,:,0], data[:,0], bounds=(l_bound, u_bound))
    model = result.x
    
    if param['source_para'] == 1:
        lnmomen = model[0]
        tstar = model[1:len(stalst)+1]
        error = model[len(stalst)+1:]
    

    if POS.upper() == "P":
        staP3lst = []
        for ista in range(len(stalst)):
            sta = stalst[ista]
            saving[sta][icase]['tstar'] = [tstar[ista]]
            saving[sta][icase]['error'] = [error[ista]]
            fitting = tstarsub.fitting(saving, sta, ORIG, POS, param['alpha'], lnmomen, 2)
            string = orid+" "+sta+" "+POS+" "+str(round(fitting, 4))
            tp.fittingfl.write(string+'\n')
            saving[sta][icase]['fitting'] = [fitting]
            if fitting >= param['misfitP']:
                staP3lst.append(sta)
        
        return staP3lst
    
    elif POS.upper() == "S":
        staS3lst = []
        for ista in range(len(stalst)):
            sta = stalst[ista]
            saving[sta][icase]['tstar'].append(tstar[ista])
            saving[sta][icase]['error'].append(error[ista])
            fitting = tstarsub.fitting(saving, sta, ORIG, POS, param['alpha'], lnmomen, 2)
            saving[sta][icase]['fitting'].append(fitting)
            string = orid+" "+sta+" "+POS+" "+str(round(fitting, 4))
            tp.fittingfl.write(string+'\n')
            if fitting >= param['misfitS']:
                staS3lst.append(sta)
        
        return staS3lst


def get_constant_correction(saving, stalst, ORIG, POS, icase, param, allsite):
    """
    #################### INVERSION FOR t* ############################
    ##    EQUATION 3 IN Stachnik, Abers, Christensen, 2004, JGR     ##
    ##      d = Gm      (Nx1) = (Nx(M+1+)) ((M+1)x1)                ##
    ##      FOR A GIVEN fc AND alpha (GRID SEARCH):                 ##
    ##      d = [ln(A1)-ln(C1)+ln(1+(f1i/fc)**2),                   ##
    ##           ln(A2)-ln(C2)+ln(1+(f2i/fc)**2),                   ##
    ##           ln(AM)-ln(CM)+ln(1+(fMi/fc)**2)]                   ##
    ##      G = [[1, -pi*f1i*f1i**(-alpha), 0, ..., 0],             ##
    ##           [1, 0, -pi*f2i*f2i**(-alpha), ..., 0],             ##
    ##           [1, 0, 0, ..., -pi*fMi*fMi**(-alpha)]]             ##
    ##      m = [[ln(Moment)],[tstar01],[tstar02],...,[tstar0M]     ##
    ##################################################################
    """

    ## icase = 3

    if POS == "P":
        bestfc = ORIG['fc']
    elif POS == "S":
        bestfc = ORIG['fc']/param['fcps']
    
    data = tstarsub.buildd(saving, stalst, ORIG, POS, icase, param, bestfc, allsite)
    G = tstarsub.buildG(saving, stalst, param['alpha'], POS, icase, param)
    
    l_bound = [-np.inf]*(2*len(stalst)+1)
    l_bound[0] = 21
    u_bound = [np.inf]*(2*len(stalst)+1)
    u_bound[0] = 40
    result = optimize.lsq_linear(G[:,:,0], data[:,0], bounds=(l_bound, u_bound))
    model = result.x
    
    if param['source_para'] == 1:
        tstar = model[1:len(stalst)+1]
        error = model[len(stalst)+1:]
    
    for ista in range(len(stalst)):
        sta = stalst[ista]
        if POS.upper() == "P":
            saving[sta][icase]['tstar'] = [tstar[ista]]
            saving[sta][icase]['error'] = [error[ista]]
            saving[sta]['updated_corr'] = [np.exp((np.log(saving[sta]['corr'][0])+saving[sta][icase]['error'][0]))]
        
        elif POS.upper() == "S":
            if len(saving[sta][icase]) == 0:
                saving[sta][icase] = {}
                saving[sta][icase]['tstar'] = [0]
                saving[sta][icase]['error'] = [1]

            saving[sta][icase]['tstar'].append(tstar[ista])
            saving[sta][icase]['error'].append(error[ista])
            if 'updated_corr' not in saving[sta]:
                saving[sta]['updated_corr'] = [1]
            
            saving[sta]['updated_corr'].append(np.exp((np.log(saving[sta]['corr'][1])+saving[sta][icase]['error'][1])))

    return saving


def inversion(orid, saving, stalst, ORIG, POS, icase, param, allsite):

    # icase = 4
    
    if POS == "P":
        bestfc = ORIG['fc']
    elif POS == "S":
        bestfc = ORIG['fc']/param['fcps']
    
    data = tstarsub.buildd(saving, stalst, ORIG, POS, icase, param, bestfc, allsite)
    G = tstarsub.buildG(saving, stalst, param['alpha'], POS, icase, param)
    try:
        Ginv = np.linalg.inv(np.dot(G[:,:,0].transpose(), G[:,:,0]))
    except:
        return ORIG, saving
    
    model, residu = optimize.nnls(G[:,:,0], data[:,0])

    if param['source_para'] == 1:
        lnmomen = model[0]             # Seismic moment
        tstar = model[1:len(stalst)+1] # t*
        momen = np.exp(lnmomen)        # Moment magnitude Mw
        Mw = float(2.0/3.0*np.log10(momen*1e7)-10.73)

        if POS == "P":
            string1 = orid+" "+POS+" "+"Mw "+str(round(Mw,2))+" "+"fc "+str(ORIG['fc'])
            tp.logfl.write(string1+"\n")
            tp.logfl.write('P wave Mw = %.2f\n' % (Mw))
            ORIG['lnmomen'] = [lnmomen] 
            ORIG['mo'] = [momen]
            ORIG['mw'] = [Mw]
        elif POS == "S":
            tp.logfl.write('S wave Mw = %.2f\n' %(Mw))
            ORIG['lnmomen'].append(lnmomen)
            ORIG['mo'].append(momen)
            ORIG['mw'].append(Mw)
    
    var = (residu**2)/(data.size-len(stalst)-2)
    ferr = open(tp.resultdir+'/%s_%serr%03d.dat' % (orid, POS.lower(), int(param['alpha']*100)),'w')
    ferr.write('%15f %7d %15f %15f\n' % (residu, data.shape[0],
                                        (residu**2)/np.sum(data[:,0]**2),
                                        (residu/np.sum(data[:,0]))))
    ferr.close()
    

    if POS.upper() == "P":
        k1 = 0
        for ista in range(len(stalst)):
            sta = stalst[ista]
            ndat = len(saving[sta][2][POS.lower()][0])
            k2 = k1 + ndat
            saving[sta][icase]['tstar'] = [tstar[ista]]
            saving[sta][icase]['misfit'] = [np.sqrt(var*(ndat-2))/ndat]    
            if param['source_para'] == 1: ## grid search for Mw, one more list for G matrix
                saving[sta][icase]['err'] = [np.sqrt(var*Ginv.diagonal()[ista+1])] ## cov(m)=cov(d)inv(G'G) FOR OVERDETERMINED PROBLEM
            else:
                saving[sta][icase]['err'] = [np.sqrt(var*Ginv.diagonal()[ista])]   ## cov(m)=cov(d)inv(G'G) FOR OVERDETERMINED PROBLEM 
            
            saving[sta][icase]['aveATTEN'] = [(1000*tstar[ista]/saving[sta]['Ptt'])]
            pfitting = tstarsub.fitting(saving, sta, ORIG, POS, param['alpha'], lnmomen, 4)
            saving[sta][icase]['fitting'] = [pfitting]

            k1 = k2
            saving[sta][icase]['Presspec'] = tstarsub.calresspec(saving[sta], POS, lnmomen, ORIG['fc'], param['alpha'])
    
    elif POS.upper() == "S":
        k1 = 0
        for ista in range(len(stalst)):
            sta = stalst[ista]
            ndat = len(saving[sta][2][POS.lower()][0])
            k2 = k1+ndat
            if len(saving[sta][icase]) == 0:
                saving[sta][icase] = {}
                saving[sta][icase]['tstar'] = [0]
                saving[sta][icase]['err'] = [1]
                saving[sta][icase]['aveATTEN'] = [0]
                saving[sta][icase]['misfit'] = [0]
                saving[sta][icase]['fitting'] = [0]
            
            saving[sta][icase]['tstar'].append(tstar[ista])
            saving[sta][icase]['misfit'].append(np.sqrt(var*(ndat-2))/ndat)
            if param['source_para'] == 1:
                saving[sta][icase]['err'].append(np.sqrt(var*Ginv.diagonal()[ista+1]))
            else:
                saving[sta][icase]['err'].append(np.sqrt(var*Ginv.diagonal()[ista]))
            saving[sta][icase]['aveATTEN'].append(1000*tstar[ista]/saving[sta]['Stt'])
            if saving[sta][icase]['aveATTEN'][0] == 0:
                saving[sta][icase]['QpQs'] = 1.75
            else:
                saving[sta][icase]['QpQs'] = saving[sta][icase]['aveATTEN'][1]/saving[sta][icase]['aveATTEN'][0]
            
            k1 = k2
            sfitting = tstarsub.fitting(saving, sta, ORIG, POS, param['alpha'], lnmomen, 4)
            saving[sta][icase]['fitting'].append(sfitting)
            saving[sta][icase]['Sresspec'] = tstarsub.calresspec(saving[sta], POS, lnmomen, ORIG['fc'], param['alpha'])

    return ORIG, saving


def output_results(orid, stalst, POS, param, ORIG, saving, ARRIV, allsite):

    if POS.upper() == "P":
        ftstar = open(tp.resultdir+'/%s_pstar%03d.dat' % (orid, int(param['alpha']*100)),'w')
        for sta in stalst:
            ftstar.write('%-12s %-4s  %.4f  %.4f  %.4f  %f  %f  %f  %f\n' %
                        (orid,
                        sta,
                        ORIG['lat'],
                        ORIG['lon'],
                        ORIG['dep'],
                        saving[sta][4]['tstar'][0],
                        saving[sta][4]['err'][0],
                        saving[sta][4]['misfit'][0],
                        saving[sta][4]['aveATTEN'][0]))
        ftstar.close()
        
        # Plot P spectrum for each station with good fitting
        if param['doplotspec']:
            for sta in stalst:
                if saving[sta][2]['good'][0]:
                    print('Plotting P spectrum of ' + sta)
                    lnM = np.log(ORIG['mo'][0])
                    tstarsub.plotspec(saving, param, ARRIV[sta], orid, 'P', lnM, ORIG['fc'], param['alpha'], 4, allsite)
        
        for sta in stalst:
            sitefl = tp.resultdir+"/site/%s_Presspec_%s.dat" %(orid, sta)
            np.savetxt(sitefl, saving[sta][4]['Presspec'], fmt='%10.4f %15.8e %6.2f')
                
    elif POS.upper() == "S":
        ftstar = open(tp.resultdir+'/%s_sstar%03d.dat' % (orid, int(param['alpha']*100)),'w')
        for sta in stalst:
            ftstar.write('%-12s %-4s  %.4f  %.4f  %.4f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n' %
                     (orid,
                      sta,
                      ORIG['lat'],
                      ORIG['lon'],
                      ORIG['dep'],
                      saving[sta][4]['tstar'][1],
                      saving[sta][4]['err'][1],
                      saving[sta][4]['misfit'][1],
                      saving[sta][4]['aveATTEN'][1],
                      saving[sta][4]['tstar'][0],
                      saving[sta][4]['err'][0],
                      saving[sta][4]['misfit'][0],
                      saving[sta][4]['aveATTEN'][0],
                      saving[sta][4]['QpQs']))
        ftstar.close()
        
        # Plot S Spectrum for each station with good fitting
        if param['doplotspec']:
            for sta in stalst:
                if saving[sta][2]['good'][1]:
                    print('Plotting S Spectrum of '+sta)
                    lnM = np.log(ORIG['mo'][1])
                    tstarsub.plotspec(saving, param, ARRIV[sta], orid, 'S', lnM, ORIG['fc'], param['alpha'], 4, allsite)
        
        for sta in stalst:
            sitefl = tp.resultdir+"/site/%s_Sresspec_%s.dat" %(orid, sta)
            np.savetxt(sitefl, saving[sta][4]['Sresspec'], fmt='%10.4f %15.8e %6.2f')

    return