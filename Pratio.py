#!/opt/antelope/5.3/bin/python
## INVERT CORNER FREQUENCY USING P WAVE SPACTRAL RATIO
## WRITTEN BY S. WEI, OCTOBER 2014, BASED ON NAKIJIMA'S FORTRAN CODES
## evt - orid:%s - mb:%f
##                 lat:%f
##                 lon:%f
##                 dep:%f
##                 stnlst:[%s]
##                 sta:%s      - freq:np.array[%f]
##                               spec:np.array[%f]

import os,glob,math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import seissign as seis

fcdir='/P/weisq/attentomo/fc_P/'
fcinvcmd='/P/weisq/attentomo/program/src/fcinvert'
specfigdir=fcdir+'specfig/'
specdir=fcdir+'spec/'
ratfigdir=fcdir+'ratiofig/'
ratdir=fcdir+'ratio/'
orlstfl=specdir+'orid.lst'
logfl=open(fcdir+'fc.detail','w')
if not os.path.isdir(ratfigdir):
    os.mkdir(ratfigdir)
if not os.path.isdir(ratdir):
    os.mkdir(ratdir)

## PARAMETERS FOR ESTIMATING RATIO
plotratio=True
shallowdep=300
##winl=11
shallowsep=30
midsep=50
deepsep=70
fcupper='15'

## SUBROUTINES
def caldist(lat1,lon1,dep1,lat2,lon2,dep2):
## Calculate the interval between two events
## USAGE: dist=caldist(lat1,lon1,dep1,lat2,lon2,dep2)
    Re=6371.0
    phi1=math.radians(90.0-lat1)
    theta1=math.radians(lon1)
    r1=Re-dep1
    x1=r1*math.cos(theta1)*math.sin(phi1)
    y1=r1*math.sin(theta1)*math.sin(phi1)
    z1=r1*math.cos(phi1)
    phi2=math.radians(90.0-lat2)
    theta2=math.radians(lon2)
    r2=Re-dep2
    x2=r2*math.cos(theta2)*math.sin(phi2)
    y2=r2*math.sin(theta2)*math.sin(phi2)
    z2=r2*math.cos(phi2)
    dist=math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    return dist


#### READ STATION INFOMATION
##stafl='/P/weisq/attentomo/input/station.lst'
##staloc={line.split()[0]:[float(line.split()[2]),float(line.split()[1])]
##     for line in open(stafl).readlines()[1:]}

## LOAD FILES FROM Pspec.py
allfreq=np.loadtxt(specdir+'allfreq')
evt={}
for line in open(orlstfl).readlines():
    evt[line.split()[0]]={}
    evt[line.split()[0]]['lat']=float(line.split()[1])
    evt[line.split()[0]]['lon']=float(line.split()[2])
    evt[line.split()[0]]['dep']=float(line.split()[3])
    evt[line.split()[0]]['mb']=float(line.split()[4])
for ordir in evt:
    stnfl=specdir+ordir+'.stn'
    stnlst=[line.rstrip() for line in open(stnfl).readlines()]
    evt[ordir]['stnlst']=stnlst
    evt[ordir]['fc']={}
    for sta in stnlst:
        freqspec=np.loadtxt(specdir+'%s_%s.spec' % (ordir,sta))
        evt[ordir][sta]={}
        evt[ordir][sta]['freq']=freqspec[:,0]
        evt[ordir][sta]['spec']=freqspec[:,1]

## INVERTING CORNER FREQUENCY FOR EACH EVENT PAIR
print('+++++++++Calculating spectra ratio for each event pair++++++++')
## LOOPING FOR EACH EVENT PAIR
for oridm in evt:     ## MASTER EVENT
    for orids in evt:     ## SLAVE EVENT
##for oridm in ['2_52102']:     ## MASTER EVENT
##    for orids in ['3_52551']:     ## SLAVE EVENT
        ## SKIP THE PAIR WITH SAME EVENT
        if oridm==orids:
            continue
        ## SKIP THE PAIR WITH MAGNITUED DIFFERENCE < 0.5
        dmg=evt[oridm]['mb']-evt[orids]['mb']
        if dmg<0.5:
            continue
        ## SKIP THE EVENT SHALLOWER THAN 50 km
        depm=evt[oridm]['dep']
        deps=evt[orids]['dep']
        if depm<50 or deps<50:
            continue
        ## SKIP THE PARI THAT ARE NOT CLUSTERED
        nocluster=True
        dist=caldist(evt[oridm]['lat'],evt[oridm]['lon'],depm,
                           evt[orids]['lat'],evt[orids]['lon'],deps)
        if (depm<shallowdep)&(deps<shallowdep)&(dist<shallowsep):
            nocluster=False
        elif (depm>shallowdep)&(deps>shallowdep)&(dist<deepsep):
            nocluster=False
        elif dist<midsep:
            nocluster=False
        if nocluster:
            continue
        print('Master: %s, Slave: %s' %(oridm,orids))
        ## INITIALIZING MATRIX FOR EACH EVENT PAIR
        ampave=np.zeros((len(allfreq))) # SUM OF log10(ampratio) FOR EACH EVT PAIR AND FREQ
        nampave=np.zeros((len(allfreq)))    # NUMBER OF AMP FOR EACH EVT PAIR AND FREQ
        ampcal=np.zeros((len(allfreq))) # SYNTHERIC LOG AVARAGE OF SOURCE SPEC RATIOS FOR ALL STATIONS
        mstn=0                          # NUMBER OF STATIONS FOR EACH EVENT PAIR
        ## AMP RATIO FOR EACH STATION PAIR FOR A PARTICULAR EVENT PAIR
        if plotratio:
            plt.figure(2)
            plt.clf()
        for stam in evt[oridm]['stnlst']:   ## STATIONS FOR MASTER EVENT
            for stas in evt[orids]['stnlst']:   ## STATIONS FOR SLAVE EVENT
##        for stam in ['A02W']:   ## STATIONS FOR MASTER EVENT
##            for stas in ['A02W']:   ## STATIONS FOR SLAVE EVENT
                if stam=='F01' or stam=='F02W':
                    continue
                if stam!=stas:
                    continue
                amp0=0.0    # SUM OF AMP RATIO FOR EACH EVENT PAIR AND STATION PAIR
                namp0=0     # NUMBER OF AMP FOR EACH EVENT PAIR AND STATION PAIR
                ampratio=np.zeros((len(allfreq)))
                freqm=evt[oridm][stam]['freq']
                freqs=evt[orids][stas]['freq']
                specm=evt[oridm][stam]['spec']
                specs=evt[orids][stas]['spec']
                for ifreq in range(len(allfreq)):
                    indm=np.nonzero(allfreq[ifreq]==freqm)[0]
                    inds=np.nonzero(allfreq[ifreq]==freqs)[0]
                    if (len(indm)*len(inds)==1):
                        if specs[inds[0]]==0:
                            ampratio[ifreq]=1000
                        else:
                            ampratio[ifreq]=specm[indm[0]]/specs[inds[0]]
                        amp0=amp0+ampratio[ifreq]
                        namp0=namp0+1
                        ampave[ifreq]=ampave[ifreq]+np.log10(ampratio[ifreq])
                        nampave[ifreq]=nampave[ifreq]+1
                ind=np.nonzero(ampratio)
##                if len(ind[0])<0:
##                    continue
##                elif len(ind[0])<winl:
##                    winl=int((len(ind[0])-1)/2)*2+1
##                ampratios=smooth(ampratio[ind])
                ampratios=ampratio[ind]
                if amp0>(2*namp0):
                    mstn=mstn+1
                if plotratio:
##                    plt.loglog(allfreq[ind],ampratio[ind],'k')
                    plt.loglog(allfreq[ind],ampratios,color='0.5')
##                    plt.text(min(allfreq[ind]),min(ampratios),stam)
                print('---Amp ratio between %s and %s' %(stam,stas))
        ## SKIP EVENT PAIR WITH LESS THAN 5 STATION PAIRS
        if mstn<5:
            continue
        for ifreq in range(len(allfreq)):
            ## ALL AVERAGE
            if nampave[ifreq]<5:
                ampave[ifreq]=0.0
            else:
                ampave[ifreq]=10**(ampave[ifreq]/nampave[ifreq])
        ind=np.nonzero(ampave)
##        ampaves=seis.smooth(ampave[ind],winl)
        ampaves=ampave[ind]
        ratiofl=ratdir+'%s.%s.ratio'%(oridm,orids)
        np.savetxt(ratiofl,np.vstack((allfreq,ampave)).transpose(),fmt='%.3f')
        cmd=fcinvcmd+' '+ratiofl+' '+fcupper
        osout=os.popen(cmd).read()
        print(osout)
        fcm=float(osout.split()[0])
        fcs=float(osout.split()[1])
        omega=float(osout.split()[2])
        if fcm==0. or fcs==0. or omega==0.:
            continue
        evt[oridm]['fc'][orids]=[fcm,dist]
        evt[orids]['fc'][oridm]=[fcs,dist]
        logfl.write('Master:%s  fc=%.2f  mb=%.1f  Slave:%s  fc=%.2f  mb=%.1f  dist=%.1f\n'
                    % (oridm,fcm,evt[oridm]['mb'],orids,fcs,evt[orids]['mb'],dist))
        ampcal=omega*((1.0+(allfreq/fcs)**2)/(1.0+(allfreq/fcm)**2))
        
        if plotratio:
            plt.loglog(allfreq[ind],ampaves,'k')
            plt.loglog(allfreq,ampcal,'r',linewidth=4)
            plt.loglog([fcm,fcm],[0.001,1000],'r--')
            plt.loglog([fcs,fcs],[0.001,1000],'r--')
            plt.text(1.5,0.5,'$mb_{%s}=%.1f:fc=%.2f$     $mb_{%s}=%.1f:fc=%.2f$' %
                     (oridm[2:],evt[oridm]['mb'],fcm,orids[2:],evt[orids]['mb'],fcs),fontsize=16)
            plt.xlim([0.5,20])
            plt.ylim([0.1,1000])
            plt.title('%s vs %s  N=%d  D=%.1f km' % (oridm,orids,mstn,dist))
            plt.savefig(ratfigdir+'%sVS%s_ratio.eps'%(oridm,orids))
logfl.close()
