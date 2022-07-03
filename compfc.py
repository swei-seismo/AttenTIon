#!/opt/antelope/5.3/bin/python
## COMPARE CORNER FREQUENCY FROM DIFFERENT METHODS
## WRITTEN BY S. WEI, OCT. 2014

import os,glob,math
import numpy as np
import matplotlib.pyplot as plt

misfcfl='/P/weisq/attentomo/tstar/bestfc.lst'
codafcfl='/P/weisq/attentomo/fc_1.5start30slogwin/fc.lst'

## IMPORT CORNER FREQUENCY BY MINIMIZING MISFIT OF t* INVERSION
misevt={}
allevt={}
for subdb in ['1','2','3','4']:
    evtlst='/P/weisq/attentomo/input/eventid_sub%s.lst' % subdb
    oridlst=[int(line.split()[0]) for line in open(evtlst)]
    maglst=[float(line.split()[1]) for line in open(evtlst)]
    for ievt in range(len(oridlst)):
        allevt[subdb+'_'+str(oridlst[ievt])]=maglst[ievt]
oridlst=[line.split()[0] for line in open(misfcfl)]
fclst=[float(line.split()[1]) for line in open(misfcfl)]
for ievt in range(len(oridlst)):
    misevt[oridlst[ievt]]=[allevt[oridlst[ievt]],fclst[ievt]]
##
##plt.figure(1)
##plt.clf()
##for dstress in [0.1,1,10,100.0]:   ## STRESS DROP IN MPa
##    mbref=np.arange(100,700)/100.0
##    Mwisc=1.54*mbref-2.54  # Mw=1.54mb-2.54, Das et al., 2011 (PREFERED)
##    mo=10**(1.5*Mwisc+9.095) ## MOMENT IN N*m
##    fcref=0.49*((dstress/mo)**(1.0/3.0))*4000*100
##    plt.plot(mbref,fcref)
##plt.legend(('$\Delta\sigma=0.1$ MPa','$\Delta\sigma=1$ MPa',
##            '$\Delta\sigma=10$ MPa','$\Delta\sigma=100$ MPa'),loc='lower left')
##plt.plot(np.array(misevt.values())[:,0],np.array(misevt.values())[:,1],
##             marker='s',linestyle='None',color='k')
##ax=plt.gca()
##ax.set_yscale('log')
##plt.xlabel('$mb_{ISC}$')
##plt.ylabel('$fc$ (Hz)')
##plt.ylim([0.05,200])
##plt.xlim([2.5,6.5])
####plt.show()
##plt.savefig('/P/weisq/attentomo/tstar/fcvsmb.eps')

## IMPORT CORNER FREQUENCY BY S-CODA RATIO
codaevt={line.split()[0]:[float(line.split()[4]),float(line.split()[1])]
         for line in open(codafcfl)}

## PLOT COMPARISON

plt.figure(2)
plt.clf()
for dstress in [0.1,1,10,100.0]:   ## STRESS DROP IN MPa
    mbref=np.arange(100,700)/100.0
    Mwisc=1.54*mbref-2.54  # Mw=1.54mb-2.54, Das et al., 2011 (PREFERED)
    mo=10**(1.5*Mwisc+9.095) ## MOMENT IN N*m
    fcref=0.49*((dstress/mo)**(1.0/3.0))*4000*100
    plt.plot(mbref,fcref)
plt.legend(('$\Delta\sigma=0.1$ MPa','$\Delta\sigma=1$ MPa',
            '$\Delta\sigma=10$ MPa','$\Delta\sigma=100$ MPa'),loc='upper right')
for evt in misevt.keys():
    try:
        codaevt[evt]
    except KeyError:
        continue
##        plt.plot(misevt[evt][0],misevt[evt][1],marker='s',linestyle='None',color='k')
    else:
        plt.plot(misevt[evt][0],misevt[evt][1],marker='s',linestyle='None',color='b')
        plt.plot(codaevt[evt][0],codaevt[evt][1],marker='s',linestyle='None',color='r')
        plt.arrow(misevt[evt][0],misevt[evt][1],codaevt[evt][0]-misevt[evt][0],
                  codaevt[evt][1]-misevt[evt][1],length_includes_head=True,
                  head_width=0.05,fc='k',ec='k',
                  head_length=10**(np.log10(codaevt[evt][1])+0.05)-10**(np.log10(codaevt[evt][1])),)
ax=plt.gca()
ax.set_yscale('log')
plt.xlabel('$mb_{ISC}$')
plt.ylabel('$fc$ (Hz)')
plt.ylim([0.1,100])
plt.xlim([2.5,6.5])
##plt.show()
plt.savefig('/P/weisq/attentomo/tstar/fccompare.eps')
