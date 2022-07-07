#!/opt/antelope/5.4/bin/python
import os
import numpy as np

s=raw_input('Please input Qnode P/S minhits ntop: ')
nodefl=s.split()[0]
pors=s.split()[1].upper()
minhits=int(s.split()[2])
ntop=int(s.split()[3])

## INPUT Q GRID
[nnx,nny,nnz]=[int(ii) for ii in open(nodefl).readlines()[0].split()]
lat=np.zeros((nnx,nny))
lon=np.zeros((nnx,nny))
nn=2
for i in range(nnx):
    for j in range(nny):
        [lat[i,j],lon[i,j]]=[float(ii) for ii in
                             open(nodefl).readlines()[nn].split()]
        nn=nn+1
##if ntop==1 or ntop==2:
dep=[float(ii) for ii in open(nodefl).readlines()[-1].split()]
##elif ntop==0:
##    nnz=nnz+1
##    dep=[12.0]
##    dep.extend([float(ii) for ii in open(nodefl).readlines()[-1].split()])

## INPUT HITS FILE
if pors=='R':
    hitsfl='hitsSP'
else:
    hitsfl='hits'+pors
hit=np.loadtxt(hitsfl)
hits=np.zeros((nnx,nny,nnz))
for i in range(nnx):
    for j in range(nny):
        for k in range(nnz):
            hits[i,j,k]=hit[i*nny+j+1-1,k]

## INPUT Q MODEL
if pors=='R' and os.path.isfile('Qinv.p'):
    qinv=np.loadtxt('Qinv.p')
    tomo=np.zeros((nnx,nny,nnz))
    nn=0
    for i in range(nnx):
        for j in range(nny):
            for k in range(nnz):
                if abs(qinv[nn]+111)<0.01 or qinv[nn]<0.01:
                    hits[i,j,k]=0
                nn=nn+1

## OUTPUT 1000/Q AT DEPTHS OF 0 AND/OR 12 KM
outfl1=open('hits'+pors+'.dep000','w')
outfl2=open('mask'+pors+'.dep000','w')
for i in range(nnx):
    for j in range(nny):
        outfl1.write('%10.4f %10.4f %10.4f\n' %
                     (lon[i,j],lat[i,j],hits[i,j,0]))
        if hits[i,j,0]<minhits:
            outfl2.write('%10.4f %10.4f NaN\n' % (lon[i,j],lat[i,j]))
        else:
            outfl2.write('%10.4f %10.4f 100\n' % (lon[i,j],lat[i,j]))
outfl1.close()
outfl2.close()
if ntop==0:
    outfl1=open('hits'+pors+'.dep012','w')
    outfl2=open('mask'+pors+'.dep012','w')
    for i in range(nnx):
        for j in range(nny):
            outfl1.write('%10.4f %10.4f %6.0f\n' %
                         (lon[i,j],lat[i,j],hits[i,j,0]))
            if hits[i,j,0]<minhits:
                outfl2.write('%10.4f %10.4f NaN\n' % (lon[i,j],lat[i,j]))
            else:
                outfl2.write('%10.4f %10.4f 100\n' % (lon[i,j],lat[i,j]))
    outfl1.close()
    outfl2.close()

## OUTPUT 1000/Q AT EACH DEPTH
for k in range(nnz):
    outfl1=open('hits%s.dep%03.0f' % (pors,dep[k]),'w')
    outfl2=open('mask%s.dep%03.0f' % (pors,dep[k]),'w')
    for i in range(nnx):
        for j in range(nny):
            outfl1.write('%10.4f %10.4f %6.0f\n' %
                         (lon[i,j],lat[i,j],hits[i,j,k]))
            if hits[i,j,k]<minhits:
                outfl2.write('%10.4f %10.4f NaN\n' % (lon[i,j],lat[i,j]))
            else:
                outfl2.write('%10.4f %10.4f 100\n' % (lon[i,j],lat[i,j]))
    outfl1.close()
    outfl2.close()
