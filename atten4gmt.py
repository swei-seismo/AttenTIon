#!/opt/antelope/5.4/bin/python
import os,glob
import numpy as np

s=raw_input('Please input resultfl Qnode P/S nnz1D ntop: ')
print(s)
qmodfl=s.split()[0]
nodefl=s.split()[1]
pors=s.split()[2].upper()
nnz1D=int(s.split()[3])
ntop=int(s.split()[4])

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

## INPUT Q MODEL
qinv=np.loadtxt(qmodfl)
tomo=np.zeros((nnx,nny,nnz))
nn=0
for i in range(nnx):
    for j in range(nny):
        for k in range(nnz):
            if abs(qinv[nn]+111)<0.01:
                qinv[nn]=-111
            elif qinv[nn]<0:
                qinv[nn]=0.00001
            if pors=='R':
                tomo[i,j,k]=qinv[nn]
            else:
                tomo[i,j,k]=qinv[nn]*1000
            nn=nn+1
qinv1D=qinv[nnx*nny*nnz:]

## INPUT PREM Q MODEL
qprem=np.loadtxt('Q.prem')
if pors=='P':
    tomoprem=1000/qprem[:,0]
elif pors=='S':
    tomoprem=1000/qprem[:,1]
elif pors=='R':
    tomoprem=qprem[:,0]/qprem[:,1]

## OUTPUT 1000/Q AT DEPTHS OF 0 AND/OR 12 KM
outfl=open('atten'+pors+'.dep000','w')
for i in range(nnx):
    for j in range(nny):
        outfl.write('%10.4f %10.4f %10.4f\n' % (lon[i,j],lat[i,j],tomoprem[0]))
outfl.close()
if ntop==0:
    outfl=open('atten'+pors+'.dep012','w')
    for i in range(nnx):
        for j in range(nny):
            outfl.write('%10.4f %10.4f %10.4f\n' % (lon[i,j],lat[i,j],tomoprem[1]))
    outfl.close()

## OUTPUT 1000/Q AT EACH DEPTH
for k in range(nnz):
    outfl=open('atten%s.dep%03.0f' % (pors,dep[k]),'w')
    for i in range(nnx):
        for j in range(nny):
            if (abs(tomo[i,j,k]+111)<0.01):
                outfl.write('%10.4f %10.4f NaN\n' %
                            (lon[i,j],lat[i,j]))
            else:
                outfl.write('%10.4f %10.4f %10.4f\n' %
                            (lon[i,j],lat[i,j],tomo[i,j,k]))
    outfl.close()

## OUTPUT Q FOR 1D MODEL
outfl=open('Q1D.result','w')
for k in range(nnz1D):
    if pors=='R':
        outfl.write('%10.4f\n' % (qinv1D[k]))
    else:
        if qinv1D[k]<=0:
            outfl.write('999999.9999\n')
        else:
            outfl.write('%10.4f\n' % (1/qinv1D[k]))
outfl.close()
