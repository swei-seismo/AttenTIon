## PROGRAM TO OUTPUT t*(P) and t*(S) with the same event-station pairs

import os
workdir='/mnt/home/swei/LauQtomo/3D1Dtomo/v2550h80top2/'
#infl='/mnt/home/swei/LauQtomo/alltstar/fcp0.5-20MPa0.27_s_tstar.dat'
infl='/mnt/home/swei/LauQtomo/tstar_fcp0.5-20MPa0.27site/s_tstar.dat'
outdir=workdir+'sameQpQs/'
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
    [lat,lon,dep,ststar,sterr,ptstar,pterr,aveQinv,QpQs]=other
    if ststar>0 and ptstar>0:
        outflp.write('%s  %.4f  %.4f  %.4f  %f  %f  %f  %.2f\n' %
                     (stn,lat,lon,dep,ptstar,pterr,pterr,aveQinv))
        outfls.write('%s  %.4f  %.4f  %.4f  %f  %f  %f  %f  %.2f  %.2f\n' %
                     (stn,lat,lon,dep,ststar,sterr,ptstar,pterr,aveQinv,QpQs))
outflp.close()
outfls.close()
