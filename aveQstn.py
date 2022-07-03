#!/Library/Frameworks/Canopy.framework/bin/python
## EXTRACT AVERAGE 1000/Q FOR EACH STATION
import numpy as np
from math import *
import matplotlib.pyplot as plt

workdir='/Users/sowei/GoogleDriveMSU/Work/Lau/Qtomo/tstar_fcp0.5-20MPa0.27site/'
stnlst='/Users/sowei/GoogleDriveMSU/Work/Lau/Qtomo/input/station.lst'

def disthead(slat,slon,flat,flon):
## Calculates distance and azimuth on sphere from starting point s
## to finishing point f
    if slon<0:
        slon=360+slon
    if flon<0:
        flon=360+flon
    slt = radians(slat)
    sln = radians(slon)
    flt = radians(flat)
    fln = radians(flon)
    delta = acos(sin(slt)*sin(flt)+cos(slt)*cos(flt)*cos(fln-sln))
    azim = atan2(sin(fln-sln)*cos(flt),
            sin(flt)*cos(slt) - cos(fln-sln)*cos(flt)*sin(slt))
    delta = degrees(delta)
    azim = degrees(azim)
    return delta,azim

## READ STATION LOCATION
stnall={}
for line in open(stnlst).readlines()[1:]:
    stnall[line.split()[0]]=[float(line.split()[1]),float(line.split()[2])]

allatten={}
## READ t* FILES AND OUTPUT
for depth in ['shallow','inter','deep']:
    allatten[depth]={}
    if depth=='shallow':
        deprange=[0,150]
    elif depth=='inter':
        deprange=[150,410]
    elif depth=='deep':
        deprange=[410,700]
#     for phase in ['p','s','r']:
    for phase in ['p']:
        aveQ={}
        aveQerr={}
        allatten[depth][phase]={}
        ## INPUT
        infl=workdir+phase+'_tstar.dat'
        if phase=='p':
            col=7
        elif phase=='s':
            col=8
        elif phase=='r':
            col=9
            infl=workdir+'s_tstar.dat'
        ievt=0
        for line in open(infl).readlines():
            dep=float(line.split()[3])
            if dep>deprange[0] and dep<deprange[1]:
                stn=line.split()[0]
                atten=float(line.split()[col])
                err=float(line.split()[5])
                tstar=float(line.split()[4])
                if tstar==0:
                    continue
                if phase=='r' and atten==1.75:
                    continue
                atterr=err/tstar*atten
                if not stn in aveQ.keys():
                    aveQ[stn]=[atten]
                    aveQerr[stn]=[atterr]
                else:
                    aveQ[stn].append(atten)
                    aveQerr[stn].append(atterr)
                ievt=ievt+1
    ##    print('%d t* measurements' % ievt)
        ## OUTPUT
        outfl=open(workdir+phase.upper()+'_stn_atten'+depth+'.dat','w')
        for stn in aveQ.keys():
    ##        print('%s has %d t* measurements' % (stn,len(aveQ[stn])))
            if len(aveQ[stn])>2:
                aveatten=np.mean(aveQ[stn])
#                 aveatten=np.median(aveQ[stn])
                stdatten=np.std(aveQ[stn])
                allatten[depth][phase][stn]=[aveatten,stdatten,len(aveQ[stn])]
        ##        aveatten=np.median(aveQ[stn])
                symsize=min(1,aveatten/stdatten/10)
                outfl.write('%4s %10.5f %10.5f %6.2f %6.2f %6.2f %d\n' %
                            (stn,stnall[stn][0],stnall[stn][1],
                             aveatten,stdatten,symsize,len(aveQ[stn])))
        outfl.close()

## PLOT AVERAGE 1/Q ALONG CROSS-SECTIONS
##plt.figure(1)
##plt.clf()
fig,axes=plt.subplots(nrows=2,ncols=3,figsize=(20,10))
plt.tight_layout(pad=4)
icol=0
for depth in ['shallow','inter','deep']:
    irow=0
    if depth=='deep':
        maxatt=10
    else:
        maxatt=25
    for xsec in ['B','A']:
        if xsec=='A':
            [rlat,rlon]=[-21.290318,183.647232]
            theta=3.5
        elif xsec=='B':
            [rlat,rlon]=[-19.730587,183.988969]
            theta=2.0
        for phase in ['p','s','r']:
            [stnlst,avelst,stdlst,distlst]=[[],[],[],[]]
            for stn in sorted(allatten[depth][phase].keys()):
                if stn.startswith(xsec):
                    stnlst.append(stn)
                    avelst.append(allatten[depth][phase][stn][0])
                    stdlst.append(allatten[depth][phase][stn][1])
                    [slon,slat]=stnall[stn]
                    (dist,az)=disthead(rlat,rlon,slat,slon)
                    dist=np.deg2rad(dist)*6367
                    if az<0:
                        dist=0-dist
                    distlst.append(dist*np.cos(np.deg2rad(theta)))
##            avelst=np.array(avelst)
##            stdlst=np.arr
            ax1=axes[irow,icol]
            if phase=='r':
                ph='Qp/Qs'
                ax2 = ax1.twinx()
                ax2.errorbar(distlst,avelst,yerr=stdlst,fmt='r-^',
                             label=ph)
                ax2.set_ylabel('Qp/Qs', color='r')
                for tl in ax2.get_yticklabels():
                    tl.set_color('r')
                ax2.set_ylim([0,5])
##                ax2.legend('Upper right')
            else:
                ph='1000/Q'+phase
                ax1.errorbar(distlst,avelst,yerr=stdlst,fmt='-o',
                             label=ph)
                ax1.set_ylabel('1000/Q')
                ax1.set_ylim([0,maxatt])
        ax1.legend(loc='upper left')
        ax2.legend(loc='upper right')
        ax1.set_xlim([-250,150])
        ax1.set_title('Stations '+xsec+': '+depth+' events')
        ax1.set_xlabel('Distance from the ridge')
        irow=irow+1
    icol=icol+1
##plt.show()
plt.savefig(workdir+'stn_xsec.eps',orientation='portrait',papertype='letter')
