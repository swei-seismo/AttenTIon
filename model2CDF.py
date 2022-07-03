#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 21:24:31 2019
Convert Q tomography models to netCDF files
@author: Shawn Wei
"""

#import math;
#import string;
import numpy as np;
import os;

def readQmod(infl,nnx,nny,nnz,dataname):
## READ INPUT Q MODEL
## USAGE: tomo=readQmod(infl,nnx,nny,nnz,modname)
    qinv=np.loadtxt(infl)
    tomo=np.zeros((nnx,nny,nnz))
    print('Read %s model: %s' % (dataname,os.path.basename(infl)))
    nn=0
    for i in range(nnx):
        for j in range(nny):
            for k in range(nnz):
#                if abs(qinv[nn]+111)<0.01:
#                    qinv[nn]=-111
#                elif qinv[nn]<0:
#                    qinv[nn]=0.00001
                if dataname=='QpQs':
                    if qinv[nn]<0:
                        tomo[i,j,k]=99999.
                    else:
                        tomo[i,j,k]=qinv[nn]
                else:
                    if qinv[nn]<=0:
                        tomo[i,j,k]=99999.
                    else:
                        tomo[i,j,k]=1./qinv[nn]
                nn=nn+1
#    qinv1D=qinv[nnx*nny*nnz:]
    return tomo


#s=input('Input model file, model name, grid file, and out dir:')
#if len(s.split())!=3:
#    exit()
#infl=s.split()[0]
#modname=s.split()[1]
#nodefl=s.split()[2]
#outdir=s.plit()p3[]

pardir=os.path.expanduser('~')+'/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQpQs/'
pardir2=os.path.expanduser('~')+'/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQps_svd/'
modname='TongaLauQ'
nodefl=pardir+'Q.grid'
outdir=os.path.expanduser('~')+'/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/netCDF/'

if not os.path.isdir(outdir):
    os.mkdir(outdir)
outfl=outdir+modname+'.metadata'
print(os.path.basename(nodefl),os.path.basename(outfl))

## Read grid file
[nnx,nny,nnz]=[int(ii) for ii in open(nodefl).readlines()[0].split()]
lat=np.zeros((nnx,nny))
lon=np.zeros((nnx,nny))
nn=2
for i in range(nnx):
    for j in range(nny):
        [lat[i,j],lon[i,j]]=[float(ii) for ii in
                             open(nodefl).readlines()[nn].split()]
        nn=nn+1
dep=[float(ii) for ii in open(nodefl).readlines()[-1].split()]

## INPUT Q MODEL
infl=pardir+'Qinv_8.0e-03_1e+04.p'
modQp=readQmod(infl,nnx,nny,nnz,'Qp')
infl=pardir+'Qinv_8.0e-03_1e+04.s'
modQs=readQmod(infl,nnx,nny,nnz,'Qk')
infl=pardir2+'Qinv_2.0_p0.9.ps'
modQpQs=readQmod(infl,nnx,nny,nnz,'QpQs')
infl=pardir2+'Qinv_2.0_p0.9.k'
modQk=readQmod(infl,nnx,nny,nnz,'Qk')

## OUTPUT CDL FILE
fout=open(outfl,'w')
fout.write('netcdf %s {\n' % modname)
fout.write('dimensions:\n')
fout.write('	depth = %d;\n' % nnz)
fout.write('	latitude = %d;\n' % nny)
fout.write('	longitude = %d;\n' % nnx)
fout.write('variables:\n')
fout.write('	float depth(depth);\n')
fout.write('            depth:long_name = \"depth below earth surface \";\n')
fout.write('            depth:units = \"kilometer\" ;\n')
fout.write('            depth:positive = \"down\" ;\n')
fout.write('	float latitude(latitude);\n')
fout.write('            latitude:long_name = \"Latitude; positive north \";\n')
fout.write('            latitude:units = \"degrees_north\" ;\n')
fout.write('            latitude:standard_name = \"latitude\" ;\n')
fout.write('	float longitude(longitude);\n')
fout.write('            longitude:long_name = \"Longitude; positive east \";\n')
fout.write('            longitude:units = \"degrees_east\" ;\n')
fout.write('            longitude:standard_name = \"Longitude\" ;\n')
fout.write('	float Qp(longitude, latitude, depth);\n')
fout.write('            Qp:long_name = \"P-wave Q\";\n')
fout.write('            Qp:missing_value = 99999.f ;\n')
fout.write('            Qp:FillValue = 99999.f ;\n')
fout.write('	float Qs(longitude, latitude, depth);\n')
fout.write('            Qs:long_name = \"S-wave Q\";\n')
fout.write('            Qs:missing_value = 99999.f ;\n')
fout.write('            Qs:FillValue = 99999.f ;\n')
fout.write('	float QpQs(longitude, latitude, depth);\n')
fout.write('            QpQs:long_name = \"Qp/Qs\";\n')
fout.write('            QpQs:missing_value = 99999.f ;\n')
fout.write('            QpQs:FillValue = 99999.f ;\n');
fout.write('	float Qk(longitude, latitude, depth);\n')
fout.write('            Qs:long_name = \"Bulk Q\";\n')
fout.write('            Qs:missing_value = 99999.f ;\n')
fout.write('            Qs:FillValue = 99999.f ;\n')
fout.write('''// global attributes:
            :title = \"%s\" ;
            :id = \"%s\" ;
            :summary = \"3D attenuation models of the Tonga-Lau mantle wedge\";
            :keywords = \"seismic, tomography, body wave attenuation, bulk attenuation, Tonga\" ;
            :Conventions = \"CF-1.0\" ;
            :Metadata_Conventions = \"Unidata Dataset Discovery v1.0\" ;
            :creator_name = \"S. Shawn Wei\" ;
            :creator_email = \"swei@msu.edu\" ;
            :repository_name = \"IRIS EMC\" ;
            :repository_url = \"http://ds.iris.edu/ds/products/emc/\" ;
            :repository_email = \"product@iris.edu\" ;
            :repository_institution = \"IRIS DMC\" ;
            :acknowledgment = \"Models were provided by S. Shawn Wei , \",
                    \"Department of Earth and Environmental Sciences\",
                    \"Michigan State University\",
                    \"St Louis, MO, 63112\" ;
            :references = \"Wei and Wiens (2019, JGR)\" ;
            :comment = \"model converted to netCDF by IRIS DMC\" ;
            :geospatial_lat_min = \"-24.0\" ;
            :geospatial_lat_max = \"-15.0\" ;
            :geospatial_lat_units = \"degrees_north\" ;
            :geospatial_lon_min = \"176\" ;
            :geospatial_lon_max = \"187\" ;
            :geospatial_lon_units = \"degrees_east\" ;
            :geospatial_vertical_min = \"25\" ;
            :geospatial_vertical_max = \"300\" ;
            :geospatial_vertical_units = \"km\" ;
            :geospatial_vertical_positive = \"down\" ;
''' % (modname,modname))

fout.write('data:\n')

# Output depth
fout.write('\n depth = ')
for i in list(range(len(dep)-1)):
    fout.write('%g, ' % (dep[i]))
fout.write('%g ;\n' % dep[-1])

# Output longitude
fout.write('\n longitude = \n')
for i in list(range(nnx)):
    for j in list(range(nny-1)):
        fout.write('%.4f, ' % (lon[i,j]))
    if i==nnx-1:
        fout.write('%.4f ;\n' % lon[i,-1])
    else:
        fout.write('%.4f,\n' % lon[i,-1])

# Output latitude
fout.write('\n latitude = \n')
for i in list(range(nnx)):
    for j in list(range(nny-1)):
        fout.write('%.4f, ' % (lat[i,j]))
    if i==nnx-1:
        fout.write('%.4f ;\n' % lat[i,-1])
    else:
        fout.write('%.4f,\n' % lat[i,-1])
# =================lat (horizontal) x lon (vertical) ==========================
# for j in list(range(lat.shape[1]))[::-1]:
#     for i in list(range(lat.shape[0]-1)):
#         fout.write('%.4f, ' % (lat[i,j]))
#     if j==0:
#         fout.write('%.4f ;\n' % lat[-1,j])
#     else:
#         fout.write('%.4f,\n' % lat[-1,j])
# =================lat (horizontal) x lon (vertical) ==========================

# Output Qp
fout.write('\n Qp = \n')
for k in list(range(nnz)):
    for i in list(range(nnx)):
        for j in list(range(nny-1)):
            fout.write('%.2f, ' % (modQp[i,j,k]))
        if k==nnz-1 and i==nnx-1:
            fout.write('%.2f ;\n' % modQp[i,-1,k])
        else:
            fout.write('%.2f,\n' % modQp[i,-1,k])
    fout.write('\n')

# Output Qs
fout.write('\n Qs = \n')
for k in list(range(nnz)):
    for i in list(range(nnx)):
        for j in list(range(nny-1)):
            fout.write('%.2f, ' % (modQs[i,j,k]))
        if k==nnz-1 and i==nnx-1:
            fout.write('%.2f ;\n' % modQs[i,-1,k])
        else:
            fout.write('%.2f,\n' % modQs[i,-1,k])
    fout.write('\n')
# =================lat (horizontal) x lon (vertical) ==========================
# for k in list(range(nnz)):
#     for j in list(range(nny))[::-1]:
#         for i in list(range(nnx-1)):
#             fout.write('%.2f, ' % (modQs[i,j,k]))
#         if k==nnz-1 and j==0:
#             fout.write('%.2f ;\n' % modQs[-1,j,k])
#         else:
#             fout.write('%.2f,\n' % modQs[-1,j,k])
#     fout.write('\n')
# =================lat (horizontal) x lon (vertical) ==========================

# Output Qp/Qs
fout.write('\n QpQs = \n')
for k in list(range(nnz)):
    for i in list(range(nnx)):
        for j in list(range(nny-1)):
            fout.write('%.2f, ' % (modQpQs[i,j,k]))
        if k==nnz-1 and i==nnx-1:
            fout.write('%.2f ;\n' % modQpQs[i,-1,k])
        else:
            fout.write('%.2f,\n' % modQpQs[i,-1,k])
    fout.write('\n')

# Output Qk
fout.write('\n Qk = \n')
for k in list(range(nnz)):
    for i in list(range(nnx)):
        for j in list(range(nny-1)):
            fout.write('%.2f, ' % (modQk[i,j,k]))
        if k==nnz-1 and i==nnx-1:
            fout.write('%.2f ;\n' % modQk[i,-1,k])
        else:
            fout.write('%.2f,\n' % modQk[i,-1,k])
    fout.write('\n')

fout.write("}");
fout.close()