#!/opt/antelope/5.4/bin/python
## PLOT MISFIT OF t* INVERSION OF
import os,glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

stnlst='/P/weisq/attentomo/input/station.lst'
evtlst='/P/weisq/attentomo/input/LauEvent.lst'
pardir='/P/weisq/attentomo/tstarPS/'
plotind=True
keydepth=70.
alphalst=['0','0.1','0.2','0.27','0.3','0.4','0.5','0.6']
dstresslst=['0.5-20','5','10','20']
fcslst=['fcp']
phaselst=['p','s']
sublst=['1','11','2','22','3','33','4','44']

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

## READ ALL t* ERROR FILES AND SAVE TO allerr, allavemis1, allavemis2
[allerr,allavemis1,allavemis2]=[{},{},{}]
for alpha in alphalst:
    allerr[alpha]={}
    allavemis1[alpha]={}
    allavemis2[alpha]={}
    for dstress in dstresslst:
        allerr[alpha][dstress]={}
        allavemis1[alpha][dstress]={}
        allavemis2[alpha][dstress]={}
        for fcs in fcslst:
            allerr[alpha][dstress][fcs]={}
            allavemis1[alpha][dstress][fcs]={}
            allavemis2[alpha][dstress][fcs]={}
            workdir=pardir+'tstar_'+fcs+dstress+'MPa'+alpha+'/'
            for phase in phaselst:
                allerr[alpha][dstress][fcs][phase]={}
                allavemis1[alpha][dstress][fcs][phase]={}
                allavemis2[alpha][dstress][fcs][phase]={}
                for subdb in sublst:
                    allavemis1[alpha][dstress][fcs][phase][subdb]=[]
                    allavemis2[alpha][dstress][fcs][phase][subdb]=[]
                    errdir=workdir+'result_sub'+subdb
                    for errfl in sorted(glob.glob(errdir+'/*'+phase+'err*.dat')):
                        evt=os.path.basename(errfl).split('_')[1]
                        errinfo=np.loadtxt(errfl)
                        allerr[alpha][dstress][fcs][phase][subdb+'_'+evt]=errinfo
                        [res,numd,norlsq,norl]=errinfo
                        dep=allevt[subdb+'_'+evt][2]
                        if dep<keydepth:
                            allavemis1[alpha][dstress][fcs][phase][subdb].append(norlsq)
                        else:
                            allavemis2[alpha][dstress][fcs][phase][subdb].append(norlsq)
            ## PLOT ERRORS VS DEPTH OR mb FOR CERTAIN PARAMETERS
            if plotind:
                fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(9,7))
                plt.tight_layout(pad=4)
                irow=0
                for phase in ['p','s']:
                    [alldep,allmb,allmis]=[[],[],[]]
                    numevt=len(allerr[alpha][dstress][fcs][phase].keys())
                    for evt in sorted(allerr[alpha][dstress][fcs][phase].keys()):
                        [dep,mb]=allevt[evt][2:4]
                        [res,numd,norlsq,norl]=allerr[alpha][dstress][fcs][phase][evt]
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
                    ax.set_title('Mean = %.1e, Std = %.1e, Nevt = %d' %
                                 (np.mean(allmis),np.std(allmis),numevt))
                    ## ERRORS VS mb
                    ax=axes[irow,1]
                    ax.scatter(allmb,allmis,c='r')
                    ax.set_xlabel('mb')
                    ax.set_ylabel('Normalized misfit of t*('+phase.upper()+')')
                    ax.set_ylim([min(allmis),max(allmis)])
                    ax.set_title('$\\alpha$ = %s, $\Delta\sigma$ = %s MPa, fcs = %s'
                                 % (alpha,dstress,fcs))
                    irow=irow+1
    ##            plt.show()
                plt.savefig(workdir+'allmisfit.eps',orientation='landscape',papertype='letter')

## PLOT TOTAL AVERAGE MISFIT OVER ALL PARAMETERS
fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(9,7))
plt.tight_layout(pad=4)
for phase in phaselst:
    if phase=='p':
        irow=0
    elif phase=='s':
        irow=1
    else:
        exit('Unknown phase: '+phase)
    for alpha in alphalst:
        for dstress in dstresslst:
            if dstress=='0.5-20':
                marker='r*'
            elif dstress=='5':
                marker='b^'
            elif dstress=='10':
                marker='ro'
            elif dstress=='20':
                marker='bs'
            else:
                exit('Unknown stress drop: '+dstress)
            for fcs in fcslst:
                [totmis1,totmis2]=[[],[]]
                for subdb in sublst:
                    totmis1.extend(allavemis1[alpha][dstress][fcs][phase][subdb])
                    totmis2.extend(allavemis2[alpha][dstress][fcs][phase][subdb])
                totmis1=np.array(totmis1)
                totmis2=np.array(totmis2)
                ## PLOT SHALLOW EVENTS
                ax=axes[irow,0]
                ax.errorbar(float(alpha),np.mean(totmis1),yerr=np.std(totmis1),
                            marker=marker[1],color=marker[0],markersize=5)
                ax.set_ylabel('Average misfit of t*('+phase.upper()+')')
##                ax.set_ylim([min(totmis2),max(totmis2)])
                ax.yaxis.get_major_formatter().set_powerlimits((0,1))
                ax.set_xlabel('$\\alpha$')
                ax.set_xlim([-0.05,0.65])
                ax.set_title('Events < %.1f km' % keydepth)
                ## PLOT DEEP EVENTS
                ax=axes[irow,1]
                ax.errorbar(float(alpha),np.mean(totmis2),yerr=np.std(totmis2),
                            marker=marker[1],color=marker[0],markersize=5)
                ax.set_ylabel('Average misfit of t*('+phase.upper()+')')
##                ax.set_ylim([min(totmis2),max(totmis2)])
                ax.yaxis.get_major_formatter().set_powerlimits((0,1))
                ax.set_xlabel('$\\alpha$')
                ax.set_xlim([-0.05,0.65])
                ax.set_title('Events > %.1f km' % keydepth)
##plt.show()
plt.savefig(pardir+'allmisfit.eps',orientation='landscape',papertype='letter')
                
#### PLOT AVERAGE MISFIT FOR EACH subdb OVER ALL PARAMETERS
fig,axes=plt.subplots(nrows=1,ncols=2,figsize=(9,7))
plt.tight_layout(pad=4)
phase='p'
for alpha in alphalst:
    for dstress in ['0.5-20']:
        for fcs in fcslst:
            for subdb in sublst:
                if '1' in subdb:
                    marker='r*'
                elif '2' in subdb:
                    marker='b^'
                elif '3' in subdb:
                    marker='ro'
                elif '4' in subdb:
                    marker='bs'
                else:
                    exit('Unknown stress drop: '+dstress)
                if len(subdb)==2:
                    fill='white'
                else:
                    fill=marker[0]
                totmis1=allavemis1[alpha][dstress][fcs][phase][subdb]
                totmis2=allavemis2[alpha][dstress][fcs][phase][subdb]
                totmis1=np.array(totmis1)
                totmis2=np.array(totmis2)
                ## PLOT SHALLOW EVENTS
                ax=axes[0]
                ax.errorbar(float(alpha),np.mean(totmis1),yerr=np.std(totmis1),
                            marker=marker[1],color=marker[0],markersize=10,
                            mfc=fill)
                ax.set_ylabel('Average misfit of t*('+phase.upper()+')')
##                ax.set_ylim([min(totmis2),max(totmis2)])
                ax.yaxis.get_major_formatter().set_powerlimits((0,1))
                ax.set_xlabel('$\\alpha$')
                ax.set_xlim([-0.05,0.65])
                ax.set_title('Events < %.1f km' % keydepth)
                ## PLOT DEEP EVENTS
                ax=axes[1]
                ax.errorbar(float(alpha),np.mean(totmis2),yerr=np.std(totmis2),
                            marker=marker[1],color=marker[0],markersize=10,
                            mfc=fill)
                ax.set_ylabel('Average misfit of t*('+phase.upper()+')')
##                ax.set_ylim([min(totmis2),max(totmis2)])
                ax.yaxis.get_major_formatter().set_powerlimits((0,1))
                ax.set_xlabel('$\\alpha$')
                ax.set_xlim([-0.05,0.65])
                ax.set_title('Events > %.1f km' % keydepth)
##plt.show()
plt.savefig(pardir+'submisfit.eps',orientation='landscape',papertype='letter')
