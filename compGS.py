#!/opt/antelope/5.4/bin/python
## PROGRAM TO COMPARE GS FILES
import os,glob
import numpy as np
import matplotlib.pyplot as plt

newdir='/P/weisq/attentomo/GS/';
olddir='/P/weisq/attentomo/backup/GS_old/';
histrange=np.hstack([0,10.**(np.arange(-3,2)),5*10.**(np.arange(-3,2))])
histrange.sort()
ratiop=[]
ratios=[]
for newpfl in sorted(glob.glob(newdir+'pgsfile*')):
    oldpfl=olddir+os.path.basename(newpfl)
    newpgs=np.loadtxt(newpfl,usecols = [0])
    oldpgs=np.loadtxt(oldpfl,usecols = [0])
    for i in range(newpgs.size):
        if oldpgs[i]==0:
            oldpgs[i]=9.e100
##            continue
        else:
            oldpgs[i]=1/oldpgs[i]
        ratiop.append((newpgs[i])/oldpgs[i])
        
for newsfl in sorted(glob.glob(newdir+'sgsfile*')):
    oldsfl=olddir+os.path.basename(newsfl)
    if os.path.getsize(newsfl)>0 and os.path.getsize(oldsfl)>0:
        newsgs=np.loadtxt(newsfl,usecols = [0])
##        print(newsfl)
        oldsgs=np.loadtxt(oldsfl,usecols = [0])
        if newsgs.size==1:
            ratios.append((newsgs)/oldsgs)
        else:
            for i in range(newsgs.size):
                if oldsgs[i]==0:
                    oldsgs[i]=9.e100
##                    continue
                else:
                    oldsgs[i]=1/oldsgs[i]
                ratios.append((newsgs[i])/oldsgs[i])

plt.figure(1)
plt.clf()
plt.subplot(2,1,1)
plt.hist(np.array(ratiop),bins=histrange)
ax=plt.gca()
ax.set_xscale('log')
plt.title('GS for P')
plt.subplot(2,1,2)
plt.hist(np.array(ratios),bins=histrange)
ax=plt.gca()
ax.set_xscale('log')
plt.xlabel('New GS/Old GS')
plt.title('GS for S')
##plt.show()
plt.savefig('/P/weisq/attentomo/GSratio.eps')

