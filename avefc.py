#!/opt/antelope/5.4/bin/python
## AVERAGING INVERTED CORNER FREQUENCY USING S-CODA WAVE RATIO
## WRITTEN BY S. WEI, JULY 2014

import os,glob,math
import numpy as np
import matplotlib.pyplot as plt

doplotave=False
doplotfcmb=True
fcdir='/P/weisq/attentomo/fc_1.5start30slogwin/'
fcinvcmd='/P/weisq/attentomo/program/src/fcinvert'
specfigdir=fcdir+'specfig/'
specdir=fcdir+'spec/'
ratfigdir=fcdir+'ratiofig/'
ratdir=fcdir+'ratio/'
orlstfl=specdir+'orid.lst'
outfl=open(fcdir+'fc.lst','w')
if not os.path.isdir(ratfigdir):
    os.mkdir(ratfigdir)
if not os.path.isdir(ratdir):
    os.mkdir(ratdir)

logfl=open(fcdir+'fc.detail','r')
evt={}
mb={}
for line in open(specdir+'orid.lst').readlines():
    mb[line.split()[0]]=float(line.split()[4])
for line in logfl.readlines():
    oridm=line.split()[0].split(':')[1]
    fcm=float(line.split()[1].split('=')[1])
    mbm=float(line.split()[2].split('=')[1])
    orids=line.split()[3].split(':')[1]
    fcs=float(line.split()[4].split('=')[1])
    mbs=float(line.split()[5].split('=')[1])
    dist=float(line.split()[6].split('=')[1])
    mb[oridm]=mbm
    mb[orids]=mbs
##    oridm=line.split()[0].split(':')[1]
##    fcm=float(line.split()[1].split('=')[1])
##    orids=line.split()[2].split(':')[1]
##    fcs=float(line.split()[3].split('=')[1])
##    dist=float(line.split()[4].split('=')[1])
    try:
        evt[oridm]
    except KeyError:
        evt[oridm]={}
    evt[oridm][orids]=[fcm,fcs,dist]
    try:
        evt[orids]
    except KeyError:
        evt[orids]={}
    evt[orids][oridm]=[fcs,fcm,dist]
logfl.close()

result=[]
resultmore=[]
for orid in evt:
    print('Corner frequency of '+orid)
    allfc=np.array(evt[orid].values())[:,0]
    alldist=np.array(evt[orid].values())[:,2]
    fcave=np.mean(allfc)
    nfc=len(allfc)
    if nfc<5:
        continue
    fcstd=np.std(allfc)
##    if fcstd==0:
##        continue
    if mb[orid]<0:
        resultmore.append([fcave,fcstd,3.0]) # No mb value
    else:
        result.append([fcave,fcstd,mb[orid]])
    if doplotave:
        plt.figure(3)
        plt.clf()
        ##plt.axhspan(fcave-fcstd,fcave+fcstd, facecolor='0.5', alpha=0.5)
        plt.scatter(alldist,allfc)
        plt.axhline(fcave,color='k',linestyle='--')
        ax=plt.gca()
        ax.set_yscale('log')
        plt.xlim([0,150])
        plt.xlabel('Inter-event distance (km)')
        plt.ylabel('Corner frequency (Hz)')
        plt.title('ORID:%s  fc=%.2f$\pm$%.2f (Hz)' % (orid[2:],fcave,fcstd))
    ##    plt.show()
        plt.savefig(ratfigdir+orid+'_fc.eps')
    outfl.write('%s   %.2f   %.2f   %d   %.1f\n' % (orid,fcave,fcstd,nfc,mb[orid]))
outfl.close()

if doplotfcmb:
    result=np.array(result)
    resultmore=np.array(resultmore)
    plt.figure(1)
    plt.clf()
    for dstress in [0.1,1,10,100.0]:   ## STRESS DROP IN MPa
        mbref=np.arange(100,700)/100.0
        Mwisc=1.54*mbref-2.54  # Mw=1.54mb-2.54, Das et al., 2011 (PREFERED)
##        Mwisc=0.85*mbref+1.03  # Mw=0.85mb+1.03, Scordilis, 2006
        mo=10**(1.5*Mwisc+9.095) ## MOMENT IN N*m
        fcref=0.49*((dstress/mo)**(1.0/3.0))*4000*100
        plt.plot(mbref,fcref)
    plt.legend(('$\Delta\sigma=0.1$ MPa','$\Delta\sigma=1$ MPa',
                '$\Delta\sigma=10$ MPa','$\Delta\sigma=100$ MPa'),loc='lower left')
    plt.errorbar(result[:,2],result[:,0],yerr=result[:,1],
                 marker='s',linestyle='None',color='k')
    plt.errorbar(resultmore[:,2],resultmore[:,0],yerr=resultmore[:,1],
                 marker='s',linestyle='None',color='r')
    ax=plt.gca()
    ax.set_yscale('log')
    plt.xlabel('$mb_{ISC}$')
    plt.ylabel('$fc$ (Hz)')
    plt.ylim([0.05,50])
    plt.xlim([2.5,6.5])
##    plt.show()
    plt.savefig(fcdir+'fcvsmb.eps')
