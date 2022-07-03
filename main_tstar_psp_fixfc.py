#!/opt/antelope/5.4/bin/python
## MAIN PROGRAM TO INVERT t*(P) and dt*(S-P) FOR ATTENUATION TOMOGRAPHY
## Written by S. Wei, May 2015, based on main_tstar_nnls.py

import os
import sys
sys.path.append(os.environ['ANTELOPE'] +'/data/python')
from antelope.datascope import *
import numpy as np
from scipy.optimize import *
from scipy.linalg import lstsq
import matplotlib.pyplot as plt
import globaldb as g
import tstarsub

## FREQUENCY DOMAIN PARAMETERS AND ALPHA
loco=0.05
hico=9.9
npol=4
##snrcrt1=[3.0,4.,2.5,2.]
snrcrt2=[2.0,3.0,2.0,2.0,1.0,1.2]   ## [SNR(P),width(P),lowSNR(P),SNR(S),width(S),lowSNR(S)]
snrcrtp1=[3.0,4.0,3.0]    ## MINIMUM [SNR,width,lowSNR] FOR FINDING BEST fc AND alpha
snrcrts1=[2.5,2.0,1.2]    ## MINIMUM [SNR,width,lowSNR] FOR FINDING BEST fc AND alpha
lincor=[0.7,0.7,0.7,1.0]    ## MINIMUM LINEAR CORRELATION COEFFICIENTS [P,S,S-P,len(S-P)]
misfitP=0.85
misfitS=0.85
##fc=np.hstack((np.arange(loco,1.01,0.05),np.arange(1.1,(hico+0.51),0.1)))
dstress=[0.5,20.0]   ## STRESS DROP IN MPa
##TIME DOMAIN PARAMETERS
WLP=5.0
WLS=8.0
pretime=50
posttime=150
prewin=[0.5,1]  ## WINDOW START BEFORE ARRIVAL
dt=0.025
beta=4000       ## SHEAR VELOCITY IN m/s FOR APPROXIMATING CORNER FREQUENCY
doplotseis=False
doplotsnr=0     ## 0 = NO SNR PLOT, 1 = ONLY S-P SNR PLOT, 2 = P & S-P SNR PLOT
doplotspec=False
doplotfcts=False

s=raw_input('Please input subdb and workdir: ')
subdb=s.split()[0]
workdir=s.split()[1]
alpha=[float(s.split()[2])]
doplotspec=bool(int(s.split()[3]))
dstress=float(s.split()[4])
##print(doplotspec)

##subdb='4'
##workdir='/P/weisq/attentomo/tstar_psp1.2_0.0'
##alpha=[0.0]

evlst='/P/weisq/attentomo/input/eventid_sub%s.lst' % subdb
oridlst=[int(line.split()[0]) for line in open(evlst).readlines()[1:]]
maglst=[float(line.split()[7]) for line in open(evlst).readlines()[1:]]
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

## IMPORT GLOBAL VARIABLES
g.dbread(subdb)

gsdir='/P/weisq/attentomo/GS'
resultdir='./result_sub'+subdb
logfile='./eventfocal%03d_sub%s.log' %(int(alpha[0]*100),subdb)
fclist=open('./bestfc.lst','a')
logfl=open(logfile,'a')
##evlst='./eventid.lst'
##oridlst=[int(line.split()[0]) for line in open(evlst)]
if not os.path.isdir(g.figdir):
    os.mkdir(g.figdir)
if not os.path.isdir(resultdir):
    os.mkdir(resultdir)

###################### LOOP OVER EACH EVENT ####################
## SAVE RESIDUL AND MISFIT OVER ALL EVENTS
allPres=0
allPmisfit=0
alltP=0
allSres=0
allSmisfit=0
alltS=0
##ORID NUMBER
for ior in range(len(oridlst)):
    orid=oridlst[ior]
##    if orid!=54501:
##        continue
    for delvar in ['Gp1','Gp2','datap1','datapf1','datap2','datapf2','resultP',
                   'Gs1','Gs2','datas1','datasf1','datas2','datasf2','resultS','findres']:
        if delvar in globals():
##            print('Deleting %s' % delvar)
            del globals()[delvar]
    logfl.write('\nWorking on ORID # %d\n' % orid)
    print('============Working on ORID # %d============' % orid)

    ## GET ORIGIN INFORMATION
    ORIG={}
    ORIG['orid']=orid
    ORIG['mb']=maglst[ior]
##    g.dbor.record=g.dbor.find('orid==%d' % (orid))
    try:
        findres=g.dbor.find('orid==%d' % (orid))
    except DbfindEnd or DbfindError:
        print('Error: No orid %d' % orid)
        continue
    else:
        g.dbor.record=findres
    
    (ORIG['time'],ORIG['lat'],ORIG['lon'],ORIG['dep']
         )=g.dbor.getv('time','lat','lon','depth')
##    (ORIG['time'],ORIG['lat'],ORIG['lon'],ORIG['dep'],ORIG['ml']
##         )=g.dbor.getv('time','lat','lon','depth','ml')

    ##DETERMINE NUMBER OF STATIONS AND RECORDS
    ##ALL P+S RECORDS LOCATED BY HUMAN
    dbsubSTA=g.dbja.subset('orid==%d && auth!~/orbassoc_l/' % orid)
    dbsubSTA=dbsubSTA.sort('sta')
    numSTA=dbsubSTA.query('dbRECORD_COUNT')
    if numSTA==0:
        print('Zero stations for ORID # %d' % orid)
        continue
    ##    exit()
    ##ALL P ARRIVELS LOCATED BY HUMAN
    dbsubP=g.dbja.subset('orid==%d && iphase=~/P/ && auth!~/orbassoc_l/' % orid)
    dbsubP=dbsubP.sort('sta')
    numP=dbsubP.query('dbRECORD_COUNT')
    ##ALL S ARRIVELS LOCATED BY HUMAN
    dbsubS=g.dbja.subset('orid==%d && iphase=~/S/ && auth!~/orbassoc_l/' % orid)
    dbsubS=dbsubS.sort('sta')
    numS=dbsubS.query('dbRECORD_COUNT')
    ##ALL S RECORDS LOCATED BY HUMAN AND ANTELOPE
    dbsubS2=g.dbja.subset('orid==%d && iphase=~/S/' % orid)
    numS2=dbsubS2.query('dbRECORD_COUNT')
    if numSTA > numS2:
        logfl.write('%d Records: %d for P and %d for S\n' % (numSTA,numP,numS2))
    elif numSTA < numS2:
        logfl.write('%d Records: %d for P and %d for S\n' % (numSTA,numP,numS2))
        logfl.write('You must have numP >= numS\n')
        continue
    ##    exit()
    else:
        logfl.write('%d Records for ORID # %d' % (numSTA,orid))

    ##LOAD FILES WITH GS (GEOMETRICAL SPREADING) AND FREE SURFACE EFFECT,
    ##    ETC FOR ZERO FREQ. AMPLITUDE
    ##print('Getting files for gemetrical spreading values...')
    pgsfile='%s/pgsfile%s_%d.txt' % (gsdir,subdb,orid)
    if not os.path.isfile(pgsfile):
        print('File %s does not exist, using mkgsfl.sh to generate' % pgsfile)
        exit()
##        os.system('/P/weisq/attentomo/program/mkgsfl.sh %s %d' % (g.dbname,orid))
    PGSTAB={'gval':[line.split()[0] for line in open(pgsfile)],'stalist':[line.split()[1] for line in open(pgsfile)]}
    sgsfile='%s/sgsfile%s_%d.txt' % (gsdir,subdb,orid)
    if not os.path.isfile(sgsfile):
        print('File %s does not exist, using mkgsfl.sh to generate' % sgsfile)
        exit()
##        os.system('/P/weisq/attentomo/program/mkgsfl.sh %s %d' % (g.dbname,orid))
    SGSTAB={'gval':[line.split()[0] for line in open(sgsfile)],'stalist':[line.split()[1] for line in open(sgsfile)]}

    ###################### LOOP OVER RECORDS IN SUBSET TABLE ####################
    dbsub=dbsubP
    ARRIV={}
    saving={}   ## saved spectra: saving[sta][icase][POS]
    staP1lst=[]     ## STATIONS USED IN FINDING BEST fc AND alpha
    staP2lst=[]     ## STATIONS USED IN t*(P) INVERSION
    staP3lst=[]     ## STATIONS USED IN t*(P) INVERSION WITHOUT BAD FITTING
    staS1lst=[]     ## STATIONS USED IN FINDING BEST fc AND alpha (NOT USED)
    staS2lst=[]     ## STATIONS USED IN t*(S) INVERSION
    staS3lst=[]     ## STATIONS USED IN t*(S) INVERSION WITHOUT BAD FITTING
    for ii in range(numP):
        dbsub.record=ii
        sta=dbsub.getv('sta')[0]
        if sta=='AFI':
            continue
##        if sta=='MSVF':
##            continue
        if sta=='F01' or sta=='F02W':
            continue
##        if sta!='B06':
##            continue
        print('Working on station %s  (%d of %d)' % (sta,ii+1,numP))

        ## DETERMINE THE CHANNEL NAMES
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
                        
        ## RETURN RAW DATA (REQUIRES P ARRIVAL, BUT NOT AN S ARRIVAL)
        ##    print("Reading seismic data from Antelope database")
        ## IF LAND STATION, dd(:,1)==TRANSVERSE && dd(:,2)==RADIAL
        ## IF OBS STATION, dd(:,1)==X && dd(:,2)==Y (x=EAST IF OBS HAS ORIENTATION IN SENSOR TABLE)
        (dd,tt,flag)=tstarsub.readseismo(pretime,dt,subdb,orid,sta,chan)
        if not flag:
            print('ERROR: Unable to read %d.%s' % (orid,sta))
            continue
        ## getdbinfo MAKES ARRIV['schan'=='EEE' OR 'EEE_00' IF NO S ARRIVAL
        (ARRIV[ii+1],Ptt,Stt,SDATA)=tstarsub.getdbinfo(sta,orid,ORIG)

        ## DETERMINE WHICH HORIZONTAL CHANNEL OF S ARRIVEL HAS LARGER AMPLITUDE
        if SDATA:
            stiim=ARRIV[ii+1]['stime']-ARRIV[ii+1]['ptime']   # S ARRIVAL IN dd SERIES
            indn=np.nonzero((tt>(stiim-1-WLS))&(tt<stiim-1))
            inds=np.nonzero((tt>stiim)&(tt<(stiim+WLS)))
            snrew=np.absolute(dd[0][inds]).max()/np.absolute(dd[0][indn]).max()
            snrns=np.absolute(dd[1][inds]).max()/np.absolute(dd[1][indn]).max()
            if snrew > snrns:
                ARRIV[ii+1]['schan']=chan[0]
            else:
                ARRIV[ii+1]['schan']=chan[1]
            
        ## WINDOW WAVEFORMS:
        ##        dd = WHOLE DATA
        ##      p_dd = P SIGNAL DATA
        ##     pn_dd = P NOISE DATA
        ##     s_dd1 = S SIGNAL DATA ON CHANNEL 1
        ##    sn_dd1 = S NOISE DATA ON CHANNEL 1
        ##    pc_dd1 = P CODA DATA ON CHANNEL 1
        ##     s_dd2 = S SIGNAL DATA ON CHANNEL 2
        ##    sn_dd2 = S NOISE DATA ON CHANNEL 2
        ##    pc_dd2 = P CODA DATA ON CHANNEL 2
        ##  PWINDATA = P SIGNAL & NOISE DATA
        ## SWINDATA1 = S SIGNAL, NOISE DATA & P CODA ON CHANNEL 1
        ## SWINDATA2 = S SIGNAL, NOISE DATA & P CODA ON CHANNEL 2
        ## IF S ARRIVEL DOES NOT EXIST, s_dd=p_dd AND sn_dd=pn_dd
        ##    print("Windowing seismic data")
        (p_dd,pn_dd,s_dd1,sn_dd1,pc_dd1,s_dd2,sn_dd2,pc_dd2,SDATA
         )=tstarsub.fixwin(dd,tt,dt,chan,ARRIV[ii+1],prewin,WLP,WLS,SDATA,
                           doplotseis,orid,sta)
        PWINDATA  = np.vstack((p_dd,pn_dd))
        SWINDATA1 = np.vstack((s_dd1,sn_dd1))
        SWINDATA1 = np.vstack((SWINDATA1,pc_dd1))
        SWINDATA2 = np.vstack((s_dd2,sn_dd2))
        SWINDATA2 = np.vstack((SWINDATA2,pc_dd2))
##        if not SDATA and max(tt)<130:
##            newln=150/dt-tt.size
##            z=np.zeros([1,int(newln)])
##            dd=np.hstack([dd,z[0]])
##            tt=np.hstack([tt,z[0]])
##        WINDATA=np.vstack((p_dd,s_dd))
##        NWINDATA=np.vstack((pn_dd,sn_dd))
##        t_dd=(np.arange(p_dd.size)+1)*dt
##        tn_dd=(np.arange(pn_dd.size)+1)*dt

        ## CALCULATE SPECTRA AND AUTO SELECTS FREQUENCY BAND ABOVE SET SNR
        ##======= 2 MEANS LOWER QUALITY DATA FOR t* INVERSION =======##
##        (spec_px,freq_px,spec_sp,freq_sp,spec,freq,n_spec,n_freq,frmn,frmx,
##         goodP2,goodS2)=tstarsub.diffspec(PWINDATA,SWINDATA1,SWINDATA2,dt,
##                        SDATA,orid,sta,snrcrtp2,snrcrts2,lincor,chan,doplotsnr)
        (spec_px,freq_px,spec,freq,n_spec,n_freq,frmn,frmx,dtstar,
         goodP2,goodS2)=tstarsub.diffspec(PWINDATA,SWINDATA1,SWINDATA2,
                SDATA,orid,sta,snrcrt2,lincor,chan,doplotsnr,alpha[0])
        if not goodP2:
##            print(sta+' No good P wave signal. Skip to next record.')
            continue
        ## SAVE SPECTRUM AND OTHER INFORMATION FOR EACH STATION
        saving[sta]={}
        saving[sta]['spec']=spec
        saving[sta]['freq']=freq
        saving[sta]['nfreq']=n_freq
        saving[sta]['nspec']=n_spec
        saving[sta]['Ptt']=Ptt
        saving[sta]['Stt']=Stt
        saving[sta][2]={}
        saving[sta][2]['good']=[goodP2,goodS2]
        saving[sta][2]['frmin']=frmn
        saving[sta][2]['frmax']=frmx
        saving[sta][2]['p']=[freq_px,spec_px]
##        saving[sta][2]['sp']=[freq_sp,spec_sp]
        ## CORRECTIONS OF GS
        correcP=float(PGSTAB['gval'][PGSTAB['stalist'].index(sta)])
        if correcP==0:
            print('Bad correction of P wave for station %s' % sta)
            continue
##        correcP=1/correcP   ## Cj
        saving[sta]['corr']=[correcP]
        staP2lst.append(sta)

##        print(sta,SDATA,goodS2,dtstar[0])
        if SDATA and goodS2 and dtstar[0]>0:
##        if SDATA and goodS2:
##            specfl=open('./'+sta+'.spec','w')
##            for ifreq in range(len(freq_sp)):
##                specfl.write('%f  %f\n' % (freq_sp[ifreq],spec_sp[ifreq]))
##            specfl.close()
            correcS=float(SGSTAB['gval'][SGSTAB['stalist'].index(sta)])
##            if correcS==0:
##                print('Bad correction of S wave for station %s' % sta)
##                continue
####            correcS=1/correcS   ## Cj
            saving[sta]['corr'].append(correcS)
            saving[sta][2]['dtstar']=dtstar
            staS2lst.append(sta)
        ##======= 2 MEANS LOWER QUALITY DATA FOR t* INVERSION =======##

        ##======= 1 MEANS HIGH QUALITY DATA FOR FINDING BEST fc AND alpha =======##        
        doplot=False
        (spec_px,freq_px,spec_sx,freq_sx,spec,freq,n_spec,n_freq,frmn,frmx,
         goodP1,goodS1)=tstarsub.dospec(PWINDATA,SWINDATA1,SWINDATA2,dt,
                        SDATA,orid,sta,snrcrtp1,snrcrts1,lincor,chan,doplot)
        if not goodP1:
##            print('No good P wave signal for finding best fc and alpha.')
            continue
        saving[sta][1]={}
        saving[sta][1]['good']=[goodP1,goodS1]
        saving[sta][1]['frmin']=frmn
        saving[sta][1]['frmax']=frmx
        saving[sta][1]['p']=[freq_px,spec_px]
        saving[sta][1]['s']=[freq_sx,spec_sx]
        staP1lst.append(sta)

        if SDATA and goodS1:
            staS1lst.append(sta)
        ##======= 1 MEANS HIGH QUALITY DATA FOR FINDING BEST fc AND alpha =======##

            
    #################### END OF LOOP OVER EACH RECORD ####################
##    logfl.write('%d records for best fc(P) and alpha, %d records for all t*(P), %d records for all t*(S)\n'
##                % (len(staP1lst),len(staP2lst),len(staS2lst)))
##    print(len(staP1lst),len(staP2lst),len(staS2lst))
    if len(staP1lst)<5:
        print('Not enough good P wave record for event %d.' % orid)
        continue
    ##    exit()


    #################### INVERSION FOR t*(P) #########################
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
        
    ## APPROXIMATE CORNER FREQUENCY RANGE BASED ON MAGNITUDE
##    ## EQUATION 4 AND 5 IN Pozgay et al, G3, 2009
##    ## CAREFUL WITH THE UNITS: fc=m/s*((10e6N/m^2)*(N*m))^(1/3)
##    ## = m/s*(10e6N/(N*m^3))^(1/3) = m/s(10^2/m) = 100/s = 100 Hz
    if ORIG['mb']<0:
        ORIG['mb']=3.0
##    ## APPROXIMATE RELATION: Mw = Ml - 0.063
##    mo=10**(1.5*(ORIG['ml']-0.063)+9.095) ## MOMENT IN N*m
    ## CONVERTING mb FROM ISC TO Mw AND MOMENT TENSOR
##    Mwisc=0.85*ORIG['mb']+1.03  # Mw=0.85mb+1.03, Scordilis, 2006
    Mwisc=1.54*ORIG['mb']-2.54  # Mw=1.54mb-2.54, Das et al., 2011 (PREFERED)
    mo=10**(1.5*Mwisc+9.095) ## MOMENT IN N*m
    bestfcp=0.49*((dstress/mo)**(1.0/3.0))*beta*100
    bestalpha=alpha[0]

    ## INVERT t*(P) WITH BEST fc AND alpha
    dataP2=tstarsub.buildd(saving,staP2lst,bestfcp,'P',2)
    GP2=tstarsub.buildG(saving,staP2lst,alpha,'P',2)
    ialco=alpha.index(bestalpha)
    GP2inv=np.linalg.inv(np.dot(GP2[:,:,ialco].transpose(),GP2[:,:,ialco]))
##    modelP,residuP,mrank,sv=lstsq(GP2[:,:,ialco],dataP2)
    modelP,residuP=nnls(GP2[:,:,ialco],dataP2[:,0])
    lnmomenP=modelP[0]      ## MOMENT
    tstarP=modelP[1:]       ## t*
    ## ESTIMATE MOMENT ERROR BASED ON ALL DATA VARIANCES
    vardatP2=residuP/np.sqrt(dataP2.shape[0]-2)
    lnmomenPerr=np.sqrt(vardatP2*GP2inv[0][0])
    ## ESTIMATE t* ERRORS BASED ON DATA VARIANCES FOR EACH t*
    estdataP=np.dot(GP2[:,:,ialco],modelP)
    staP3lst=[]
    k1=0
    for ista in range(len(staP2lst)):
        sta=staP2lst[ista]
        ndat=len(saving[sta][2]['p'][0])
        k2=k1+ndat
        dat=dataP2[k1:k2]
        est=estdataP[k1:k2]
##        saving[sta][2]['est']=[est]
##        saving[sta][2]['dat']=[dat]
        var=(np.linalg.norm(dat-est)**2)/(ndat-2)    ## POSTERIOR VARIANCE USED AS PRIOR VARIANCE
        saving[sta][2]['tstar']=[tstarP[ista]]
        saving[sta][2]['misfit']=[np.sqrt(var*(ndat-2))/ndat]
        saving[sta][2]['err']=[np.sqrt(var*GP2inv.diagonal()[ista+1])] ## cov(m)=cov(d)inv(G'G) FOR OVERDETERMINED PROBLEM
        saving[sta][2]['aveATTEN']=[(1000*tstarP[ista]/saving[sta]['Ptt'])]
        if saving[sta][2]['good'][0]:
            pfitting=tstarsub.fitting(saving[sta],sta,orid,'P',lnmomenP,
                              bestfcp,bestalpha,2)
            saving[sta][2]['fitting']=[pfitting]
        else:
            saving[sta][2]['fitting']=[1000]
##        if np.sqrt(var)<1:
##        if saving[sta][2]['err'][0]<0.03:
##        if saving[sta][2]['misfit'][0]<misfitP:
        if saving[sta][2]['fitting'][0]>misfitP:
            staP3lst.append(sta)
        k1=k2

###### 3 MEANS INVERTING AGAIN WITHOUT misfit>1 ######
    if len(staP3lst)<5:
        print('Not enough good P wave record for event %d.' % orid)
        continue        
    dataP3=tstarsub.buildd(saving,staP3lst,bestfcp,'P',3)
    GP3=tstarsub.buildG(saving,staP3lst,alpha,'P',3)
    ialco=alpha.index(bestalpha)
    GP3inv=np.linalg.inv(np.dot(GP3[:,:,ialco].transpose(),GP3[:,:,ialco]))
##    modelP,residuP,mrank,sv=lstsq(GP3[:,:,ialco],dataP3)
    modelP,residuP=nnls(GP3[:,:,ialco],dataP3[:,0])
    lnmomenP=modelP[0]      ## MOMENT
    tstarP=modelP[1:]       ## t*
    ## SAVE RESIDUL AND MISFIT OVER ALL
    allPres=allPres+(residuP**2)/dataP2.shape[0]
    allPmisfit=allPmisfit+residuP/dataP2.shape[0]
    alltP=alltP+1
    ferr=open(resultdir+'/%s_%d_perr%03d.dat' % (subdb,orid,int(alpha[0]*100)),'w')
    ferr.write('%15f %7d %15f %15f\n' % (residuP,dataP2.shape[0],
                                     ((residuP**2)/dataP2.shape[0]),
                                     (residuP/(dataP2.shape[0]))))
    ferr.close()
    ## ESTIMATE MOMENT ERROR BASED ON ALL DATA VARIANCES
    vardatP3=residuP/np.sqrt(dataP3.shape[0]-2)
    lnmomenPerr=np.sqrt(vardatP3*GP3inv[0][0])
    ## ESTIMATE t* ERRORS BASED ON DATA VARIANCES FOR EACH t*
    estdataP=np.dot(GP3[:,:,ialco],modelP)
    k1=0
    for ista in range(len(staP3lst)):
        sta=staP3lst[ista]
        ndat=len(saving[sta][2]['p'][0])
        k2=k1+ndat
        dat=dataP3[k1:k2]
        est=estdataP[k1:k2]
        var=(np.linalg.norm(dat-est)**2)/(ndat-2)    ## POSTERIOR VARIANCE USED AS PRIOR VARIANCE
        saving[sta][3]={}
##        saving[sta][3]['est']=[est]
##        saving[sta][3]['dat']=[dat]
        saving[sta][3]['tstar']=[tstarP[ista]]
        saving[sta][3]['misfit']=[np.sqrt(var*(ndat-2))/ndat]
        saving[sta][3]['err']=[np.sqrt(var*GP3inv.diagonal()[ista+1])] ## cov(m)=cov(d)inv(G'G) FOR OVERDETERMINED PROBLEM
        saving[sta][3]['aveATTEN']=[(1000*tstarP[ista]/saving[sta]['Ptt'])]
        k1=k2
###### 3 MEANS INVERTING AGAIN WITHOUT misfit>1 ######
    
    ## MOMENT MAGNITUDE Mw
    momenP=np.exp(lnmomenP)
    Mw=float(2.0/3.0*np.log10(momenP*1e7)-10.73)
    logfl.write('Mw = %.2f, Mwisc = %.2f\n' % (Mw,Mwisc))

    ## OUTPUT RESULTS FOR TOMOGRAPHY
    ftstar=open(resultdir+'/%s_%d_pstar%03d.dat' % (subdb,orid,int(alpha[0]*100)),'w')
    for icntp in range(len(staP3lst)):
        sta=staP3lst[icntp]
        ftstar.write('%s  %.4f  %.4f  %.4f  %f  %f  %f  %.2f\n' % (sta,
                ORIG['lat'],ORIG['lon'],ORIG['dep'],saving[sta][3]['tstar'][0],
                saving[sta][3]['err'][0],saving[sta][3]['misfit'][0],
                saving[sta][3]['aveATTEN'][0]))
    ftstar.close()

    ## PLOT P SPECTRUM FOR EACH STATION
    if doplotspec:
        for sta in staP3lst:
            if saving[sta][2]['good'][0]:
                print('Plotting P spectrum of ' + sta)
                tstarsub.plotspec(saving[sta],sta,orid,'P',lnmomenP,
                                  bestfcp,bestalpha,3)
    

    ## OUTPUT dt*(S-P) AND t*(S)
    if len(staS2lst)<1:
        print('Not enough good S wave record for event %d.' % orid)
        continue
    ftstar=open(resultdir+'/%s_%d_spstar%03d.dat' % (subdb,orid,int(alpha[0]*100)),'w')
    for icnts in range(len(staS2lst)):
        sta=staS2lst[icnts]
        if sta in staP3lst:
            dtstar=saving[sta][2]['dtstar'][0]
            sterr=saving[sta][2]['dtstar'][1]
            ptstar=saving[sta][3]['tstar'][0]
            pterr=saving[sta][3]['err'][0]
            ststar=ptstar+dtstar
            patten=saving[sta][3]['aveATTEN'][0]
            satten=1000*ststar/saving[sta]['Stt']
            allSres=allSres+saving[sta][2]['dtstar'][2]
            allSmisfit=allSmisfit+saving[sta][2]['dtstar'][3]
            alltS=alltS+1
            if patten==0:
                QpQs=1.75
            else:
                QpQs=satten/patten
            ftstar.write('%s  %.4f  %.4f  %.4f  %f  %f  %f  %f  %f  %.2f\n' %
                         (sta,ORIG['lat'],ORIG['lon'],ORIG['dep'],ststar,
                          sterr,ptstar,patten,satten,QpQs))
    ftstar.close()


## AVERAGE RESIDUAL AND MISFIT FOR P AND S
alltP=max(alltP,1)
alltS=max(alltS,1)
logfl.write('%f %f %d %f %f %d\n' % (allPres,allPmisfit,alltP,allSres,allSmisfit,alltS))
logfl.write('%f %f %d %f %f %d\n' % (allPres/alltP,allPmisfit/alltP,alltP,allSres/alltS,allSmisfit/alltS,alltS))

logfl.close()
fclist.close()
