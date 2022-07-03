#!/opt/antelope/5.4/bin/python
## READ Qp AND Qs MODEL INDEPENDENTLY AND CALCULATE Qp/Qs
import os,shutil
import numpy as np
from scipy.interpolate import *

pardir='/P/weisq/attentomo/3D1Dtomo/'
Qpdir=pardir+'v26h30/fcp0.5-20MPa0.27/'
Qsdir=pardir+'v70h80/fcp0.5-20MPa0.27/'
Qpmodfl=open(Qpdir+'Qinv_3.0e-04_1e+04.p','r')
Qsmodfl=Qsdir+'Qinv_2.8e-02_1e+04.s'
Qpnodfl=open(Qpdir+'Q.grid','r')
Qsnodfl=open(Qsdir+'Q.grid','r')
Qs1D=4

QpQsdir=pardir+'v70h80/indQpQs/'
if not os.path.isdir(QpQsdir):
    os.mkdir(QpQsdir)
shutil.copyfile(Qsdir+'Q.grid',QpQsdir+'Q.grid')
shutil.copyfile(Qsdir+'Q.prem',QpQsdir+'Q.prem')
shutil.copyfile(Qsdir+'hitsS',QpQsdir+'hitsSP')
QpQsmodfl=QpQsdir+'QpQs.sp'

def lonlat2xy(lat0,lon0,beta,lat1,lon1):
##  CONVERT lat, lon to CARTISIAN COORDINATES RELATIVE TO lat0, lon0
    re=6377
    if lon1<0:
        lon1 = 360.+lon1
    if lon0<0:
        lon0 = 360.+lon0
    xx = np.deg2rad(lon1-lon0)*re
    yy = np.deg2rad(lat1-lat0)*re
    x1 = (xx-yy*np.tan(beta))*np.cos(beta)
    y1 = x1*np.tan(beta)+yy/np.cos(beta)
    return x1,y1

## IMPORT Qs NODES
[nx,ny,nz]=[int(i) for i in Qsnodfl.readline().split()]
[reflat,reflon,beta]=[float(i) for i in Qsnodfl.readline().split()]
Qx=np.zeros(nx)
Qy=np.zeros(ny)
for i in range(nx):
    for j in range(ny):
        [lat,lon]=[float(irow) for irow in Qsnodfl.readline().split()]
        (Qx[i],Qy[j])=lonlat2xy(reflat,reflon,beta,lat,lon)
Qz=np.array([float(k) for k in Qsnodfl.readline().split()])
Qsnode=np.empty([nx*ny*nz,3])
inode=0
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            Qsnode[inode]=[Qx[i],Qy[j],Qz[k]]
            inode=inode+1

## IMPORT Qs MODEL
Qsinv=np.loadtxt(Qsmodfl)
Qsinv=Qsinv[0:-Qs1D]
##Qsinv=np.empty((nx,ny,nz))
##for i in range(nx):
##    for j in range(ny):
##        for k in range(nz):
##            Qsinv[i,j,k]=float(Qsmodfl.readline())

## IMPORT Qp NODES
[nx,ny,nz]=[int(i) for i in Qpnodfl.readline().split()]
[tmplat,tmplon,beta]=[float(i) for i in Qpnodfl.readline().split()]
Qx=np.zeros(nx)
Qy=np.zeros(ny)
for i in range(nx):
    for j in range(ny):
        [lat,lon]=[float(irow) for irow in Qpnodfl.readline().split()]
        (Qx[i],Qy[j])=lonlat2xy(reflat,reflon,beta,lat,lon)
Qz=np.array([float(k) for k in Qpnodfl.readline().split()])
Qpnode=(Qx,Qy,Qz)

## IMPORT Qp MODEL
Qpinv=np.empty((nx,ny,nz))
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            Qpinv[i,j,k]=float(Qpmodfl.readline())

## INTERPOLATE Qp TO Qs GRID
##print('%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' %
##      (max(Qpnode[0]),min(Qpnode[0]),max(Qpnode[1]),
##       min(Qpnode[1]),max(Qpnode[2]),min(Qpnode[2])))
##print('%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' %
##      (max(Qsnode[:,0]),min(Qsnode[:,0]),max(Qsnode[:,1]),
##       min(Qsnode[:,1]),max(Qsnode[:,2]),min(Qsnode[:,2])))
Qpinvnew=interpn(Qpnode,Qpinv,Qsnode,bounds_error=False,fill_value=None)

## CALCULATE Qp/Qs
QpQs=Qsinv/Qpinvnew
print(max(QpQs),np.mean(QpQs),np.median(QpQs),min(QpQs))
QpQs=np.append(QpQs,np.zeros(3))
np.savetxt(QpQsmodfl,QpQs,fmt='%12.7e')
