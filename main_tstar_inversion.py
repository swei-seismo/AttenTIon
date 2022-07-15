#!/usr/bin/env python3
# -*- coding:utf-8 -*-
## MAIN PROGRAM TO INVERT t* FOR ATTENUATION TOMOGRAPHY
## Written by S. Wei, Nov. 2013
## Edited by Yurong Zhang, Jun. 2022

import os
from importlib_metadata import files
import tstar_parameters as tp
import tstar_load
import tstar_inversion_function as tf
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size, rank = comm.Get_size(), comm.Get_rank()

param = tp.set_parameters()
oridlst = tp.working_dir()

## problem: data preprocessing? deal with sac files with picked arrival times

num = int(len(oridlst)/size)
bgn = rank*num
num_last = len(oridlst) - (size-1)*num
NUM = (rank<size-1 and num or num_last)

for ior in np.arange(bgn,bgn+NUM,1):

    orid=oridlst[ior]
    tp.logfl.write('\nWorking on ORID # %s\n' % orid)
    print('============Working on # %s (%d of %d)============' % (orid,ior+1,len(oridlst)))
    
    resultfl=tp.resultdir+'/%s_pstar%03d.dat' % (orid,int(param["alpha"]*100))
    if os.path.isfile(resultfl):
        print('Skip %s, already exists' % orid)
        continue
    
    ##Load event and station information
    (stalst, ORIG, ARRIV) = tstar_load.loaddata(param, orid, tp.sacdir+"/%s"%orid)
    if stalst == 0:
        print('Zero stations for # %s' % orid)
        continue
    os.chdir(tp.workdir)

    ##Load constants with geometric spreading and free surface effects
    (PGS, SGS) = tstar_load.loadGS(orid, tp.gsdir)
    ##Loop over records
    (staP1lst, staP2lst, staP3lst, staS1lst, staS2lst, staS3lst, saving, ARRIV) = \
        tf.Loop_Record(orid, stalst, param, ARRIV, PGS, SGS)
    ##Find the best fc if grid searching
    if param['source_para'] == 1:
        if len(staP1lst)<5:
            print('Not enough good P wave record for event %s.' % orid)
            continue
        ##grid search for the best fc
        (ORIG, flag) = tf.bestfc(orid, saving, staP1lst, ORIG, 'P', 1, param)
        if flag == 0:
            continue

    ##2nd inversion: INVERT t*(P) WITH BEST fc
    (ORIG, saving) = tf.inversion(orid, saving, staP2lst, ORIG, 'P', 2, param)
    for sta in staP2lst:
        if saving[sta][2]['fitting'][0]>param['misfitP']:
            staP3lst.append(sta)
    if len(staP3lst)<5:
        print('Not enough good P wave record for event %s.' % orid)
        continue  
    ##3rd inversion: NVERTING AGAIN WITHOUT fitting < 0.85
    (ORIG, saving) = tf.inversion(orid, saving, staP3lst, ORIG, 'P', 3, param)

    ##output the results
    tf.output_results(orid, staP3lst, param, ORIG, saving)


tp.logfl.close()
tp.fclist.close()


