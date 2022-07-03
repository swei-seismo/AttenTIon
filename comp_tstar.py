#%%
# Compare t* measurements from different inversions

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Helvetica','font.size': 36})

infl='/Users/swei/GoogleDriveMSU/Work/Lau/Qtomo/newQ/tstar/p_tstar.dat'
stnfl='/Users/swei/GoogleDriveMSU/Work/Lau/Qtomo/input/station.lst'
maxdep=300  # Max depth of intermediate-depth earthquakes
keyevt=[-17.7675,  -178.0745,   620.9924]  # 1_51844
keyevt2=[-22.1087,  -179.3094,   598.4839]  # 4_54521

# Import station list
stnlst=pd.read_csv(stnfl,sep='\s+',header=0,
                   names=['stnname','stnlon','stnlat','stnele','num','orientation'])

# Import tstar results
tstardata=pd.read_csv(infl,sep='\s+',header=0,
                      names=['stnname','evtlat','evtlon','evtdep','newt','newt_err','newaveatten','newaveatten_err','oldt','oldt_err','oldaveatten','oldaveatten_err'])

# # Add station info to tstar list
# tstardata['stnlat']=np.nan
# tstardata['stnlon']=np.nan
# tstardata['stnele']=np.nan
# for ii in range(len(tstardata)):
#     tstardata.stnlat[ii]=stnlst.stnlat[stnlst.stnname==tstardata.stnname[ii]]
#     tstardata.stnlon[ii]=stnlst.stnlon[stnlst.stnname==tstardata.stnname[ii]]
#     tstardata.stnele[ii]=stnlst.stnele[stnlst.stnname==tstardata.stnname[ii]]

# Plot
fig,axs=plt.subplots(nrows=2,ncols=2,figsize=(30,20))
fig.tight_layout(pad=4)
# Plot t* comparison
ax=axs[0,0]
ax.set(title='t*(P)',aspect=1,
           xlim=(0,max([max(tstardata.oldt),max(tstardata.newt)])),
           ylim=(0,max([max(tstardata.oldt),max(tstardata.newt)])),
           xlabel='WW2018 t* (s)',ylabel='New t* (s)')
# axs[0].scatter(tstardata.oldt,tstardata.newt,s=50,
#                facecolors='none',edgecolors='b')
ax.scatter(tstardata.oldt[tstardata.evtdep>maxdep],
               tstardata.newt[tstardata.evtdep>maxdep],
               s=50,facecolors='none',edgecolors='b',alpha=0.5,
               label='Event depth > %d km' % maxdep)
ax.scatter(tstardata.oldt[tstardata.evtdep<=maxdep],
               tstardata.newt[tstardata.evtdep<=maxdep],
               s=200,facecolors='none',edgecolors='r',alpha=0.5,
               label='Event depth <= %d km' % maxdep)
# axs[0].scatter(tstardata.oldt[tstardata.stnele>0],
#                tstardata.newt[tstardata.stnele>0],
#                s=50,facecolors='none',edgecolors='b')
# axs[0].scatter(tstardata.oldt[tstardata.stnele<=0],
#                tstardata.newt[tstardata.stnele<=0],
#                s=200,facecolors='none',edgecolors='r',alpha=0.1)
ax.plot([0,100],[0,100],'k--',linewidth=2)
ax.legend(loc='lower right', fontsize=22)
# Plot histgram
ax=axs[1,0]
histrange=(-100,100)
ax.set(xlabel='t* change (%)',ylabel='Count')
ax.hist((tstardata.newt[tstardata.evtdep>maxdep]-
         tstardata.oldt[tstardata.evtdep>maxdep])/tstardata.oldt[tstardata.evtdep>maxdep]*100,
        bins=20,range=histrange,color='b',alpha=0.5,label='Event depth > %d km' % maxdep)
ax.hist((tstardata.newt[tstardata.evtdep<=maxdep]-
         tstardata.oldt[tstardata.evtdep<=maxdep])/tstardata.oldt[tstardata.evtdep<=maxdep]*100,
        bins=20,range=histrange,color='r',alpha=0.5,label='Event depth <= %d km' % maxdep)
ax.legend(loc='upper left', fontsize=22)
# Plot path-average attenuation comparison
ax=axs[0,1]
ax.set(title='Path-average P-wave attenuation',aspect=1,
           xlim=(0,max([max(tstardata.oldaveatten),max(tstardata.newaveatten)])),
           ylim=(0,max([max(tstardata.oldaveatten),max(tstardata.newaveatten)])),
           xlabel='WW2018 path-average 1/Q',ylabel='New path-average 1/Q')
# axs[1].scatter(tstardata.oldaveatten,tstardata.newaveatten,s=50,
#                facecolors='none',edgecolors='b')
ax.scatter(tstardata.oldaveatten[tstardata.evtdep>maxdep],
               tstardata.newaveatten[tstardata.evtdep>maxdep],
               s=50,facecolors='none',edgecolors='b',alpha=0.5,
               label='Event depth > %d km' % maxdep)
ax.scatter(tstardata.oldaveatten[tstardata.evtdep<=maxdep],
               tstardata.newaveatten[tstardata.evtdep<=maxdep],
               s=200,facecolors='none',edgecolors='r',alpha=0.5,
               label='Event depth <= %d km' % maxdep)
# axs[1].scatter(tstardata.oldaveatten[tstardata.stnele>0],
#                tstardata.newaveatten[tstardata.stnele>0],
#                s=50,facecolors='none',edgecolors='b')
# axs[1].scatter(tstardata.oldaveatten[tstardata.stnele<=0],
#                tstardata.newaveatten[tstardata.stnele<=0],
#                s=200,facecolors='none',edgecolors='r',alpha=0.1)
ax.plot([0,100],[0,100],'k--',linewidth=2)
ax.legend(loc='lower right', fontsize=22)
# Plot histgram
ax=axs[1,1]
ax.set(xlabel='Path-average 1/Q change (%)',ylabel='Count')
ax.hist((tstardata.newaveatten[tstardata.evtdep>maxdep]-
         tstardata.oldaveatten[tstardata.evtdep>maxdep])/tstardata.oldaveatten[tstardata.evtdep>maxdep]*100,
        bins=20,range=histrange,color='b',alpha=0.5,label='Event depth > %d km' % maxdep)
ax.hist((tstardata.newaveatten[tstardata.evtdep<=maxdep]-
         tstardata.oldaveatten[tstardata.evtdep<=maxdep])/tstardata.oldaveatten[tstardata.evtdep<=maxdep]*100,
        bins=20,range=histrange,color='r',alpha=0.5,label='Event depth <= %d km' % maxdep)
ax.legend(loc='upper left', fontsize=22)
# # Plot 1_51844
# keyind=np.nonzero([np.logical_and(np.logical_and(tstardata.evtlon==keyevt[1],tstardata.evtlat==keyevt[0]),tstardata.evtdep==keyevt[2])])[1]
# axs[0,0].set(title='t*(P) of 1_51844',aspect=1,
#            xlim=(0,max([max(tstardata.oldt),max(tstardata.newt)])),
#            ylim=(0,max([max(tstardata.oldt),max(tstardata.newt)])),
#            xlabel='WW2018 t* (s)',ylabel='New t* (s)')
# axs[0,0].scatter(tstardata.oldt[keyind],
#                tstardata.newt[keyind],
#                s=200,facecolors='none',edgecolors='b')
# axs[0,0].plot([0,100],[0,100],'k--',linewidth=2)
# # Plot 4_54521
# keyind=np.nonzero([np.logical_and(np.logical_and(tstardata.evtlon==keyevt2[1],tstardata.evtlat==keyevt2[0]),tstardata.evtdep==keyevt2[2])])[1]
# axs[0,1].set(title='t*(P) of 4_54521',aspect=1,
#            xlim=(0,max([max(tstardata.oldt),max(tstardata.newt)])),
#            ylim=(0,max([max(tstardata.oldt),max(tstardata.newt)])),
#            xlabel='WW2018 t* (s)',ylabel='New t* (s)')
# axs[0,1].scatter(tstardata.oldt[keyind],
#                tstardata.newt[keyind],
#                s=200,facecolors='none',edgecolors='b')
# axs[0,1].plot([0,100],[0,100],'k--',linewidth=2)
# # # Plot path-average attenuation comparison
# # axs[1,1].set(title='Path-average P-wave attenuation of 1_51844',aspect=1,
# #            xlim=(0,max([max(tstardata.oldaveatten),max(tstardata.newaveatten)])),
# #            ylim=(0,max([max(tstardata.oldaveatten),max(tstardata.newaveatten)])),xlabel='WW2018 path-average 1/Q',ylabel='New path-average 1/Q')
# # axs[1,1].scatter(tstardata.oldaveatten[keyind],
# #                tstardata.newaveatten[keyind],
# #                s=200,facecolors='none',edgecolors='b')
# # axs[1,1].plot([0,100],[0,100],'k--',linewidth=2)
# fig.show()
# %%
