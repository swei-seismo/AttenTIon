## PROGRAM TO OUTPUT t*(P) and t*(S) with the same event-station pairs

import os
#workdir='/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/v125h155top2/'
#if not os.path.isdir(workdir):
#    os.mkdir(workdir)
infl='/Users/swei/GoogleDriveMSU/Work/Lau/Qtomo/alltstar_nosite/fcp0.5-20MPa0.27_s_tstar.dat'
#outdir=workdir+'sameQpQs_fcp/'
outdir='/Users/swei/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQpQs/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
outflp=open(outdir+'p_tstar.dat','w')
outfls=open(outdir+'s_tstar.dat','w')

## READ ORIGINAL t*(S) FILE AND OUTPUT
for line in open(infl).readlines():
    stn=line.split()[0]
    other=line.split()[1:]
    for i in range(len(other)):
        other[i]=float(other[i])
    [lat,lon,dep,ststar,sterr,ptstar,pterr,aveQsinv,QpQs]=other
    if ststar>0 and ptstar>0:
        outflp.write('%s  %.4f  %.4f  %.4f  %f  %f  %f  %.2f\n' %
                     (stn,lat,lon,dep,ptstar,pterr,pterr,aveQsinv/QpQs))
        outfls.write('%s  %.4f  %.4f  %.4f  %f  %f  %f  %f  %.2f  %.2f\n' %
                     (stn,lat,lon,dep,ststar,sterr,ptstar,pterr,aveQsinv,QpQs))
outflp.close()
outfls.close()
