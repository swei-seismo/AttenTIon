#!/opt/antelope/5.4/bin/python
## COMPARE MISFIT OF t* INVERSION FOR EACH EVENT
import os,glob
import numpy as np
import matplotlib.pyplot as plt

alpha='0.6'
dstress='0.5-20'
fcs='fcp'
workdir='/P/weisq/attentomo/tstarPS/tstar_'+fcs+dstress+'MPa'+alpha+'/'
stnlst='/P/weisq/attentomo/input/station.lst'
evtlst='/P/weisq/attentomo/input/LauEvent.lst'

## READ EVENT INFOMATION
allevt={}
for line in open(evtlst).readlines()[1:]:
    [lauid,laulat,laulon,laudep,date,hms,mb,iscid]=line.split()
    [subdb,orid]=lauid.split('_')
    laulat=float(laulat)
    laulon=float(laulon)
    if laulon<0:
        laulon=laulon+360 
    laudep=float(laudep)
    mb=float(mb)
    allevt[lauid]=[laulat,laulon,laudep,mb]

## READ ALL t* ERROR FILES
allerr={}
for phase in ['p','s']:
    allerr[phase]={}
    for subdb in ['1','11','2','22','3','33','4','44']:
        errdir=workdir+'result_sub'+subdb
        for errfl in sorted(glob.glob(errdir+'/*'+phase+'err*.dat')):
            evt=os.path.basename(errfl).split('_')[1]
            errinfo=np.loadtxt(errfl)
            allerr[phase][subdb+'_'+evt]=errinfo

## PLOT ERRORS VS DEPTH OR mb
fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(9,7))
plt.tight_layout(pad=4)
irow=0
for phase in ['p','s']:
    [alldep,allmb,allmis]=[[],[],[]]
    numevt=len(allerr[phase].keys())
    for evt in sorted(allerr[phase].keys()):
        [dep,mb]=allevt[evt][2:4]
        [res,numd,norlsq,norl]=allerr[phase][evt]
        alldep.append(dep)
        allmb.append(mb)
        allmis.append(norlsq)
    alldep=np.array(alldep)
    allmb=np.array(allmb)
    allmis=np.array(allmis)
    ## ERRORS VS DEPTH
    ax=axes[irow,0]
    ax.scatter(alldep,allmis,c='b')
    ax.set_xlabel('Depth (km)')
    ax.set_ylabel('Normalized misfit of t*('+phase.upper()+')')
    ax.set_ylim([min(allmis),max(allmis)])
##    ax.set_title('%d events for %s wave' % (numevt,phase.upper()))
    ax.set_title('Mean = %f, Std = %f' % (np.mean(allmis),np.std(allmis)))
    ## ERRORS VS mb
    ax=axes[irow,1]
    ax.scatter(allmb,allmis,c='r')
    ax.set_xlabel('mb')
    ax.set_ylabel('Normalized misfit of t*('+phase.upper()+')')
    ax.set_ylim([min(allmis),max(allmis)])
##    ax.set_title('%d events for %s wave' % (numevt,phase.upper()))
    ax.set_title('$\\alpha$ = %s, $\Delta\sigma$ = %s MPa, fcs = %s'
                 % (alpha,dstress,fcs))
    irow=irow+1
plt.show()
##plt.savefig(workdir+'allmisfit.eps',orientation='landscape',papertype='letter')
