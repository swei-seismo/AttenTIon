import os,glob
import numpy as np

s=raw_input('Please input resultfl Qnode P/S nnz1D ntop: ')
print(s)
qmodfl=s.split()[0]
nodefl=s.split()[1]
pors=s.split()[2]
nnz1D=int(s.split()[3])
ntop=int(s.split()[4])

if len(pors)==1:
    pors=pors.upper()

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
            if pors=='R' or pors=='QmQk':
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
elif pors=='QmQk':
    tomoprem=qprem[:,0]/qprem[:,1]
    tomoprem[0]=600./57823.
    tomoprem[1]=600./57823.
elif pors=='Qk':
    tomoprem=qprem[:,0]/qprem[:,1]
    tomoprem[0]=1000./57823.
    tomoprem[1]=1000./57823.

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
#            outfl.write('%10.4f %10.4f %10.4f\n' % (lon[i,j],lat[i,j],tomoprem[1]))
            if (abs(tomo[i,j,k]+111)<0.01) and pors=='R':
                outfl.write('%10.4f %10.4f NaN\n' %
                            (lon[i,j],lat[i,j]))
            elif (abs(tomo[i,j,k]+111)<0.01) and pors=='QmQk':
                outfl.write('%10.4f %10.4f 0.01\n' %
                            (lon[i,j],lat[i,j]))
            elif (abs(tomo[i,j,k]+111000)<0.01) and pors=='Qk':
                outfl.write('%10.4f %10.4f 0.0173\n' %
                            (lon[i,j],lat[i,j]))
            else:
                outfl.write('%10.4f %10.4f %10.4f\n' %
                            (lon[i,j],lat[i,j],tomo[i,j,k]))
    outfl.close()

## OUTPUT 1000/Q AT EACH DEPTH
for k in range(nnz):
    outfl=open('atten%s.dep%03.0f' % (pors,dep[k]),'w')
    for i in range(nnx):
        for j in range(nny):
            if (abs(tomo[i,j,k]+111)<0.01) and pors=='R':
                outfl.write('%10.4f %10.4f NaN\n' %
                            (lon[i,j],lat[i,j]))
            elif (abs(tomo[i,j,k]+111)<0.01) and pors=='QmQk':
                outfl.write('%10.4f %10.4f 0.01\n' %
                            (lon[i,j],lat[i,j]))
            elif (abs(tomo[i,j,k]+111000)<0.01) and pors=='Qk':
                outfl.write('%10.4f %10.4f 0.0173\n' %
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

## READ AND OUTPUT MODEL RESOLUTION AND VARIANCE IF EXIST
modresfl=qmodfl+'_varm'
if os.path.isfile(modresfl) and pors=='R':
    qinvvar=np.loadtxt(modresfl)
    tomovar=np.zeros((nnx,nny,nnz))
    tomores=np.zeros((nnx,nny,nnz))
    nn=0
    for i in range(nnx):
        for j in range(nny):
            for k in range(nnz):
                tomores[i,j,k]=qinvvar[nn,0]
                tomovar[i,j,k]=qinvvar[nn,1]
                nn=nn+1
    qinvvar1D=qinvvar[nnx*nny*nnz:]
    for k in range(nnz):
        outfl=open('atten%svar.dep%03.0f' % (pors,dep[k]),'w')
        outfl2=open('atten%sres.dep%03.0f' % (pors,dep[k]),'w')
        for i in range(nnx):
            for j in range(nny):
                outfl.write('%10.4f %10.4f %10.4f\n' %
                            (lon[i,j],lat[i,j],tomovar[i,j,k]))
                outfl2.write('%10.4f %10.4f %10.4f\n' %
                            (lon[i,j],lat[i,j],tomores[i,j,k]))
        outfl.close()
elif os.path.isfile(modresfl) and (pors=='P' or pors=='S' or pors=='Qk'):
    qinvvar=np.loadtxt(modresfl)
    tomovar=np.zeros((nnx,nny,nnz))
    nn=0
    for i in range(nnx):
        for j in range(nny):
            for k in range(nnz):
                tomovar[i,j,k]=qinvvar[nn]*1000
                nn=nn+1
    qinvvar1D=qinvvar[nnx*nny*nnz:]
    for k in range(nnz):
        outfl=open('atten%sstd.dep%03.0f' % (pors,dep[k]),'w')
        for i in range(nnx):
            for j in range(nny):
                outfl.write('%10.4f %10.4f %10.4f\n' %
                            (lon[i,j],lat[i,j],tomovar[i,j,k]))

        outfl.close()
