#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:36:59 2019
Convert Q tomography models to GeoCSV files
@author: Shawn Wei
"""

#import math;
#import string;
import numpy as np
import os,datetime

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
                elif dataname=='QpQsstd':
                    tomo[i,j,k]=qinv[nn,1]
                elif 'inv' in dataname:
                    if qinv[nn]<=0:
                        tomo[i,j,k]=0.
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


# pardir=os.path.expanduser('~')+'/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQpQs/'
# pardir2=os.path.expanduser('~')+'/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQps_svd/'
# modname='TongaLau.Q.2019'
pardir=os.path.expanduser('~')+'/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v25h30top2_Qp/fcp0.5-20MPa0.27/'
# pardir2=os.path.expanduser('~')+'/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQps_svd/'
modname='TongaLau.Qp.2018'
nodefl=pardir+'Q.grid'
outdir=os.path.expanduser('~')+'/GoogleDriveMSU/Work/Lau/Qtomo/IRISEMC/'

if not os.path.isdir(outdir):
    os.mkdir(outdir)
outfl=outdir+modname+'.csv'
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
infl=pardir+'Qinv_7.0e-04_1e+06.p'
modQp=readQmod(infl,nnx,nny,nnz,'Qp')
# infl=pardir+'Qinv_8.0e-03_1e+04.s'
# modQs=readQmod(infl,nnx,nny,nnz,'Qs')
# infl=pardir2+'Qinv_2.0_p0.9.ps'
# modQpQs=readQmod(infl,nnx,nny,nnz,'QpQs')
# infl=pardir2+'Qinv_2.0_p0.9.k'
# modQk=readQmod(infl,nnx,nny,nnz,'Qk')
# infl=pardir+'Qinv_8.0e-03_1e+04.p_varm'
# modQpstd=readQmod(infl,nnx,nny,nnz,'invQpstd')
# infl=pardir+'Qinv_8.0e-03_1e+04.s_varm'
# modQsstd=readQmod(infl,nnx,nny,nnz,'invQsstd')
# infl=pardir2+'Qinv_2.0_p0.9.ps_varm'
# modQpQsstd=readQmod(infl,nnx,nny,nnz,'QpQsstd')
# infl=pardir2+'Qinv_2.0_p0.9.k_varm'
# modQkstd=readQmod(infl,nnx,nny,nnz,'invQkstd')

## OUTPUT CDL FILE
now = datetime.datetime.utcnow()
fout=open(outfl,'w')
fout.write('''# dataset: GeoCSV2.0
# created: %s UTC
# delimiter: |
# global_title: 3D P-wave attenuation models of the Tonga-Lau mantle wedge
# global_id: %s
# global_summary: High-resolution 3D model of Qp of the Tonga-Lau mantle wedge
# global_reference: Wei and Wiens (2018, EPSL)
# global_references: http://ds.iris.edu/ds/products/emc-earthmodels/
# global_keywords: seismic, tomography, body wave attenuation, P-wave attenuation, Tonga
# global_Conventions: CF-1.0
# global_Metadata_Conventions: Unidata Dataset Discovery v1.0
# global_creator_name: IRIS EMC
# global_creator_url: http://ds.iris.edu/ds/products/emc/
# global_creator_email: product@iris.edu
# global_institution: IRIS DMC
# global_acknowledgment: Models were provided by S. Shawn Wei, Department of Earth and Environmental Sciences, Michigan State University
# global_comment: model converted to GeoCSV by S. Wei
# global_geospatial_lat_min: -24.0
# global_geospatial_lat_max: -15.0
# global_geospatial_lat_units: degrees_north
# global_geospatial_lon_min: 176
# global_geospatial_lon_max: 187
# global_geospatial_lon_units: degrees_east
# global_geospatial_vertical_min: 25
# global_geospatial_vertical_max: 300
# global_geospatial_vertical_units: km
# global_geospatial_vertical_positive: down
# latitude_column: latitude
# latitude_long_name: Latitude; positive north
# latitude_units: degrees_north
# latitude_standard_name: latitude
# longitude_column: longitude
# longitude_long_name: Longitude; positive east
# longitude_units: degrees_east
# longitude_standard_name: longitude
# depth_column: depth
# depth_long_name: depth below earth surface
# depth_units: km
# depth_positive: down
# Qp_column: Qp
# Qp_long_name: P-wave Q
# Qp_display_name: Qp
# Qp_units: count
# Qp_missing_value: 99999.0
# Qp__FillValue: 99999.0
longitude | latitude | depth | Qp
''' % (now.strftime('%Y-%m-%d %H:%M:%S'),modname))

# Output models
for k in list(range(nnz)):
    for i in list(range(nnx)):
        for j in list(range(nny)):
#             fout.write('%9.4f | %9.4f | %5.1f | %9.2f | %9.2f | %5.2f | %9.2f | %9.4f | %9.4f | %5.2f | %9.4f\n'
#                        % (lon[i,j],lat[i,j],dep[k],
#                           modQp[i,j,k],modQs[i,j,k],modQpQs[i,j,k],modQk[i,j,k],
#                           modQpstd[i,j,k],modQsstd[i,j,k],modQpQsstd[i,j,k],modQkstd[i,j,k]))
            fout.write('%9.4f | %9.4f | %5.1f | %9.2f\n'
                       % (lon[i,j],lat[i,j],dep[k],modQp[i,j,k]))
fout.close()