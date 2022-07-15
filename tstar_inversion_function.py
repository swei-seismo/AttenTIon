#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##Main function for t* inversion 
import numpy as np
from scipy.optimize import *
from scipy.linalg import lstsq
from matplotlib import pyplot as plt
import tstar_parameters as tp
import tstar_load 
import tstarsub

def Loop_Record(orid, stalst, param, ARRIV, PGS, SGS):
    """
    Loop over each record, process the seismic information

    :param staP1lst:

    :param staP2lst:

    """
    saving={}   ## saved spectra: saving[sta][icase][POS]
    staP1lst=[]     ## STATIONS USED IN FINDING BEST fc AND alpha
    staP2lst=[]     ## STATIONS USED IN t*(P) INVERSION
    staP3lst=[]     ## STATIONS USED IN t*(P) INVERSION WITHOUT BAD FITTING
    staS1lst=[]     ## STATIONS USED IN FINDING BEST fc AND alpha (NOT USED)
    staS2lst=[]     ## STATIONS USED IN t*(S) INVERSION
    staS3lst=[]     ## STATIONS USED IN t*(S) INVERSION WITHOUT BAD FITTING
    for i, sta in enumerate(stalst):
        if sta=='AFI':
            continue
        if sta=='F01' or sta=='F02W':
            continue 
        ##yurong: question: are they bad stations?
        print('Working on station %s  (%d of %d)' % (sta,i+1,len(stalst)))
        chan = tstar_load.channels(sta)

        ## RETURN RAW DATA (REQUIRES P ARRIVAL, BUT NOT AN S ARRIVAL)
        ##    print("Reading seismic data from Antelope database")
        ## IF LAND STATION, dd(:,1)==TRANSVERSE && dd(:,2)==RADIAL
        ## IF OBS STATION, dd(:,1)==X && dd(:,2)==Y (x=EAST IF OBS HAS ORIENTATION IN SENSOR TABLE)
        (dd,tt,flag)=tstarsub.readseismo(param['pretime'],param['dt'],orid,sta,chan)
        if not flag:
            print('ERROR: Unable to read %s.%s' % (orid,sta))
            continue

        ## DETERMINE WHICH HORIZONTAL CHANNEL OF S ARRIVEL HAS LARGER AMPLITUDE
        if ARRIV[sta]['SDATA']:
            stiim=ARRIV[sta]['T1']-ARRIV[sta]['T0']   # S ARRIVAL IN dd SERIES
            indn=np.nonzero((tt>(stiim-1-param['WLS']))&(tt<stiim-1))
            ##yurong: questions: why -1 in NS direction?
            inds=np.nonzero((tt>stiim)&(tt<(stiim+param['WLS'])))
            snrew=np.absolute(dd[0][inds]).max()/np.absolute(dd[0][indn]).max()
            snrns=np.absolute(dd[1][inds]).max()/np.absolute(dd[1][indn]).max()
            if snrew > snrns:
                ARRIV.setdefault(sta,{})['schan'] = chan[0]
            else:
                ARRIV.setdefault(sta,{})['schan'] = chan[1]

        (p_dd,pn_dd)=tstarsub.fixwin(dd,tt,param,chan,ARRIV[sta],orid,sta)
        ## yurong: work on S wave
        PWINDATA  = np.vstack((p_dd,pn_dd))

        ## CALCULATE SPECTRA AND AUTO SELECTS FREQUENCY BAND ABOVE SET SNR
        ##======= 2 MEANS LOWER QUALITY DATA FOR t* INVERSION =======##
        (goodP2, goodS2, spec_px, freq_px, pspec, pn_spec, pfreq, pn_freq, frmin, frmax) \
            = tstarsub.dospec(PWINDATA, orid, sta, param, chan, 2)
        if not goodP2:
            print('No good P wave signal. Skip to next record.')
            continue

        ## SAVE SPECTRUM AND OTHER INFORMATION FOR EACH STATION
        saving[sta]={}
        saving[sta]['spec'], saving[sta]['nspec'] = [pspec], [pn_spec]
        saving[sta]['freq'], saving[sta]['nfreq'] = [pfreq], [pn_freq]
        saving[sta][1], saving[sta][2], saving[sta][3] = {}, {}, {}

        saving[sta]['Ptt'] = ARRIV[sta]['T0']
        saving[sta][2]['good']=[goodP2,goodS2]
        saving[sta][2]['frmin'], saving[sta][2]['frmax'] = frmin, frmax
        saving[sta][2]['p']=[freq_px,spec_px]
        # saving[sta][2]['s']=[freq_sx,spec_sx]

        ## CORRECTIONS OF GS
        correcP=float(PGS['gval'][PGS['stalist'].index(sta)])
        if correcP==0:
            print('Bad correction of P wave for station %s' % sta)
            continue
        saving[sta]['corr']=[correcP]
        staP2lst.append(sta)

        if ARRIV[sta]['SDATA'] and goodS2:
            correcS=float(SGS['gval'][SGS['stalist'].index(sta)])
            if correcS==0:
                print('Bad correction of S wave for station %s' % sta)
                continue
            saving[sta]['corr'].append(correcS)
            staS2lst.append(sta)
        ##======= 2 MEANS LOWER QUALITY DATA FOR t* INVERSION =======##

        if param['source_para'] == 1:   ## search for fc
            ##======= 1 MEANS HIGH QUALITY DATA FOR FINDING BEST fc AND alpha =======##        
            doplot=False
            # (spec_px,freq_px,spec_sx,freq_sx,spec,freq,n_spec,n_freq,frmn,frmx,
            # goodP1,goodS1)=tstarsub.dospec(PWINDATA,SWINDATA1,SWINDATA2,dt,
            #                 SDATA,orid,sta,snrcrtp1,snrcrts1,lincor,chan,doplot)
            (goodP1, goodS1, spec_px, freq_px, pspec, pn_spec, pfreq, pn_freq, frmin, frmax) \
                 = tstarsub.dospec(PWINDATA, orid, sta, param, chan, 1)
            if not goodP1:
                # print('No good P wave signal for finding best fc and alpha.')
                continue
            saving[sta][1]['good']=[goodP1,goodS1]
            saving[sta][1]['frmin']=frmin
            saving[sta][1]['frmax']=frmax
            saving[sta][1]['p']=[freq_px,spec_px]
            # saving[sta][1]['s']=[freq_sx,spec_sx]
            staP1lst.append(sta)

            if ARRIV[sta]['SDATA'] and goodS1 and goodS2:
                staS1lst.append(sta)
            ##======= 1 MEANS HIGH QUALITY DATA FOR FINDING BEST fc AND alpha =======##

    return staP1lst, staP2lst, staP3lst, staS1lst, staS2lst, staS3lst,\
        saving, ARRIV

def inversion(orid, saving, stalst, ORIG, POS, icase, param):
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
    idata = (icase==3 and 2 or icase)
    
    data = tstarsub.buildd(saving,stalst,ORIG,POS,idata,param['source_para'],ORIG['fc'])
    G = tstarsub.buildG(saving,stalst,param['alpha'],POS,idata,param['source_para'])
    Ginv=np.linalg.inv(np.dot(G[:,:,0].transpose(),G[:,:,0]))
    model,residu=nnls(G[:,:,0],data[:,0])
    if param['source_para'] == 1:
        lnmomen=model[0]      ## MOMENT
        tstar=model[1:]       ## t*
        ## MOMENT MAGNITUDE Mw
        momen=np.exp(lnmomen)
        Mw=float(2.0/3.0*np.log10(momen*1e7)-10.73)
        ORIG['mo'] = momen
        ORIG['mw'] = Mw
        # tp.logfl.write('Mw = %.2f, Mwisc = %.2f\n' % (Mw,Mwisc))
    elif param['source_para'] == 2:
        tstar=model           ## t*
        lnmomen=np.log(ORIG['mo'])
    elif param['source_para'] == 3:
        tstar=model           ## t*
        # lnmomen=np.log(ORIG['mo'])
    if icase == 2:
        ferr=open(tp.resultdir+'/%s_perr%03d.dat' % (orid,int(param['alpha']*100)),'w')
        ferr.write('%15f %7d %15f %15f\n' % (residu,data.shape[0],
                                        (residu**2)/np.sum(data[:,0]**2),
                                        (residu/np.sum(data[:,0]))))
        ferr.close()

    ## ESTIMATE MOMENT ERROR BASED ON ALL DATA VARIANCES
    vardat=residu/np.sqrt(data.shape[0]-2)
    # lnmomenPerr=np.sqrt(vardatP2*GP2inv[0][0])
    ## ESTIMATE t* ERRORS BASED ON DATA VARIANCES FOR EACH t*
    estdata=np.dot(G[:,:,0],model)
    var = (residu**2)/(data.size-len(stalst)-2)
    
    k1=0
    for ista in range(len(stalst)):
        sta = stalst[ista]
        ndat=len(saving[sta][idata][POS.lower()][0])
        k2=k1+ndat
        dat=data[k1:k2]
        est=estdata[k1:k2]
        # var=(np.linalg.norm(dat-est)**2)/(ndat-2)    ## POSTERIOR VARIANCE USED AS PRIOR VARIANCE
        
        saving[sta][icase]['tstar']=[tstar[ista]]
        saving[sta][icase]['misfit']=[np.sqrt(var*(ndat-2))/ndat]    
        if param['source_para'] == 1: ## grid search for Mw, one more list for G matrix
            saving[sta][icase]['err']=[np.sqrt(var*Ginv.diagonal()[ista+1])] ## cov(m)=cov(d)inv(G'G) FOR OVERDETERMINED PROBLEM
        else:
            saving[sta][icase]['err']=[np.sqrt(var*Ginv.diagonal()[ista])] ## cov(m)=cov(d)inv(G'G) FOR OVERDETERMINED PROBLEM        saving[sta][icase]['aveATTEN']=[(1000*tstar[ista]/saving[sta]['Ptt'])]
        saving[sta][icase]['aveATTEN']=[(1000*tstar[ista]/saving[sta]['Ptt'])]
        if icase == 2:
            ## Measure how synthetic curve fit the observed data
            if saving[sta][icase]['good'][0]:
                pfitting=tstarsub.fitting(saving,sta,ORIG,POS,param['alpha'],lnmomen,2)
                saving[sta][icase]['fitting']=[pfitting]
            else:
                saving[sta][icase]['fitting']=[1000]
        elif icase == 3:
            lnM = np.log(ORIG['mo'])
            saving[sta][icase]['resspec']=[tstarsub.calresspec(saving[sta],POS,lnM,
                                                ORIG['fc'],param['alpha'])]
        k1=k2

    return ORIG, saving

def bestfc(orid, saving, stalst, ORIG, POS, icase, param):
    """ APPROXIMATE CORNER FREQUENCY RANGE BASED ON MAGNITUDE
    EQUATION 4 AND 5 IN Pozgay et al, G3, 2009
    CAREFUL WITH THE UNITS: fc=m/s*((10e6N/m^2)*(N*m))^(1/3)
    = m/s*(10e6N/(N*m^3))^(1/3) = m/s(10^2/m) = 100/s = 100 Hz
    """
    if ORIG['mb']<=0:
        ORIG['mb']=3
    Mwisc=1.54*ORIG['mb']-2.54  # Mw=1.54mb-2.54, Das et al., 2011 (PREFERED)
    mo=10**(1.5*Mwisc+9.095) ## MOMENT IN N*m
    fclow=0.49*((param['dstress'][0]/mo)**(1.0/3.0))*param['beta']*100
    fchigh=0.49*((param['dstress'][1]/mo)**(1.0/3.0))*param['beta']*100
    if fclow<1 and fchigh<=1.1:
        fc=np.arange(fclow,fchigh,0.02)
    elif fclow<1 and fchigh>1.1:
        fc=np.hstack((np.arange(fclow,1.09,0.02),np.arange(1.1,fchigh,0.1)))
    else:
        fc=np.arange(fclow,fchigh,0.1)
    if max(fc)<fchigh:
        fc=np.hstack((fc,fchigh))
    if fc.shape[0]<5:
        tp.logfl.write('Too narrow frequency band for finding fc')
        return ORIG, 0
    tp.logfl.write('fc for 0.5 MPa: %.2f,fc for 20 MPa: %.2f\n' % (fclow,fchigh))
    tp.logfl.write('fc(P) will range from %.2f to %.2f for P\n' % (min(fc),max(fc)))
    
    ## BUILD G MATRIX TO FIND BEST fc AND alpha
    G = tstarsub.buildG(saving,stalst,param['alpha'],POS,icase,param['source_para'])
    Ginv=np.linalg.inv(np.dot(G[:,:,0].transpose(),G[:,:,0]))
    tsfc=np.zeros((len(fc),len(stalst)))
    for ifc in range(len(fc)):
        data = tstarsub.buildd(saving,stalst,ORIG,POS,icase,param['source_para'],fc[ifc])
        model,residu = nnls(G[:,:,0],data[:,0])
        lnmomen = model[0]      ## MOMENT
        tstar = model[1:]       ## t*
        L2P = residu/np.sum(data[:,0])
        vardat = L2P/(data.shape[0]-len(stalst)-1)
        lnmomen_err = np.sqrt(vardat*Ginv[0][0])
        tstar_err = np.sqrt(vardat*Ginv.diagonal()[1:])
        tsfc[ifc,:]=tstar
        try:
            result
        except NameError:
            result=np.array([[fc[ifc],lnmomen,L2P,vardat,lnmomen_err]])
        else:
            result=np.vstack((result,np.array([[fc[ifc],lnmomen,L2P,vardat,lnmomen_err]])))
    L2Pall=result[:,2].tolist()
    bestresult=result[L2Pall.index(min(L2Pall))]
    bestfc=float(bestresult[0])
    ORIG['fc'] = bestfc
    tp.logfl.write('Best fc(P) = %.2f Hz\n' % (bestfc))
    if max(fc)==bestfc:
        tp.logfl.write('Warning: best fc is upper limit of fc\n')
    if bestfc==min(fc):
        tp.logfl.write('Warning: best fc is lower limit of fc\n')
    tp.fclist.write('%s   %.2f  %.1f\n' % (orid,bestfc,ORIG['mb']))

    ## PLOTTING L2P VS CORNER FREQUENCY
    if param['doplotfcall']:
        fig=plt.figure(10)
        fig.clf()
        fig.subplots_adjust(wspace=0.3,hspace=0.3)
        ax1=fig.add_subplot(1,2,1)
        ax1.plot(result[:,0],L2Pall,'b*-')
        ax1.plot(bestfc,min(L2Pall),'r^',ms=10)
        ax1.set_xlabel('Corner Frequency (Hz)')
        ax1.set_ylabel('L2 Norm')
        ax1.set_title('ORID = %s' % (orid))
        ax2=fig.add_subplot(1,2,2)
        ax2.plot(result[:,0],np.log10(np.exp(result[:,1])*1e7),'b*-')
        ax2.plot(bestfc,np.log10(np.exp(bestresult[1])*1e7),'r^',ms=10)
        ax2.set_xlabel('Corner Frequency (Hz)')
        ax2.set_ylabel('log10(moment)')
        ax2.set_title('ORID = %s' % (orid))
        fig.savefig(tp.figdir5+'/%s_fcP.pdf' % (orid))
        
    ## PLOT t*(P) VS CORNER FREQUENCY FOR EACH STATION
    if param['doplotfcts']:
        for ista in range(len(stalst)):
            plt.figure(20)
            plt.clf()
            tspert=(tsfc[:,ista]-np.mean(tsfc[:,ista]))/np.mean(tsfc[:,ista])*100
            plt.plot(fc,tspert,'b*-')
            plt.plot(bestfc,tspert[fc==bestfc],'r^',ms=10)
            plt.title('ORID = %s at %s' % (orid,stalst[ista]))
            plt.xlabel('fc')
            plt.ylabel('t* perturbation (%)')
            plt.savefig(tp.figdir4+'/%s_%s_fcts.pdf' % (orid,stalst[ista]))
 

    return ORIG, 1

def output_results(orid, staP3lst, param, ORIG, saving):
    ## OUTPUT P RESIDUAL SPECTRA FOR SITE EFFECTS
    for sta in staP3lst:
        sitefl=tp.resultdir+'/%s_Presspec_%s.dat' % (orid,sta)
        np.savetxt(sitefl,saving[sta][3]['resspec'][0], fmt='%10.4f  %15.8e  %6.2f')

    ## OUTPUT RESULTS FOR TOMOGRAPHY
    ftstar=open(tp.resultdir+'/%s_pstar%03d.dat' % (orid,int(param['alpha']*100)),'w')
    for sta in staP3lst:
        ftstar.write('%s  %.4f  %.4f  %.4f  %f  %f  %f  %.2f\n' %
                     (sta,
                      ORIG['lat'],
                      ORIG['lon'],
                      ORIG['dep'],
                      saving[sta][3]['tstar'][0],
                      saving[sta][3]['err'][0],
                      saving[sta][3]['misfit'][0],
                      saving[sta][3]['aveATTEN'][0]))
    ftstar.close()

    ## PLOT P SPECTRUM FOR EACH STATION
    if param['doplotspec']:
        for sta in staP3lst:
            if saving[sta][2]['good'][0]:
                print('Plotting P spectrum of ' + sta)
                lnM = np.log(ORIG['mo'])
                print(lnM,ORIG['fc'],ORIG['mo'])
                tstarsub.plotspec(saving[sta],sta,orid,'P',lnM,
                                  ORIG['fc'],param['alpha'],3)

    return

