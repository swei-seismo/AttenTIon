## PROGRAM TO CREATE GEOMETRIC SPREADING FILES FOR P AND S
## Written by S. Wei, MARCH 2015

import os
import numpy as np
import momcalc
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

main_path = os.getcwd()
workdir = main_path+"/GS"
namedir = main_path+"/Eventname"
catalog = main_path+"/AACSE_catalog.dat"
sacdir = main_path+"/correctedSeismograms"
if not os.path.exists(workdir):
        os.makedirs(workdir, exist_ok=True)

oridlst = []
for name in open(namedir).readlines():
    name = name.strip('\n')
    if name is not None:
        oridlst.append(name)

num = int(len(oridlst)/size)
bgn = rank*num
num_last = len(oridlst)-(size-1)*num
NUM = (rank<size-1 and num or num_last)

for ior in np.arange(bgn,bgn+NUM,1):
    orid = oridlst[ior]
    print(f"Rank {rank}: {orid}")
    command = "cat %s | awk '/%s/ {print}'" %(catalog, orid)
    line = os.popen(command).read().strip()
    depth = float(line.split()[3])

    outpfln = workdir+'/pgsfile_%s.txt' % (orid)
    outsfln = workdir+'/sgsfile_%s.txt' % (orid)
    outpfile = open(outpfln,'w')
    outsfile = open(outsfln,'w')

    dbev = {}
    dstime = {}
    
    command = "saclst KSTNM t1 gcarc f %s/%s/*Z.sac" %(sacdir,orid)
    output = os.popen(command).read().strip().split("\n")
    for line in output:
        line = line.strip()
        sta = line.split()[1]
        stime = line.split()[2]
        gcarc = line.split()[3]
        dbev[sta] = float(gcarc)
        dstime[sta] = float(stime)

    dbev = sorted(dbev.items(), key=lambda item:item[1])
    num = len(dbev)
    
    inimomp = np.array([1.0e-12, 5.0e-13, 1.0e-13])
    inimoms = np.array([1.0e-12, 5.0e-13, 1.0e-13])
    for i in range(num):
        sta = dbev[i][0]
        delta = dbev[i][1]
        depth = max(depth,0.00)  
        for phs in ['P', 'S']:
            moment = momcalc.correct(delta, depth, phs)
            if phs == 'P':
                if moment != 0:
                    inimomp[0] = inimomp[1]
                    inimomp[1] = inimomp[2]
                    inimomp[2] = moment
                if moment == 0:
                    moment = np.mean(inimomp[np.nonzero(inimomp)])
            elif phs == 'S' and dstime[sta] != -12345:
                if moment != 0:
                    inimoms[0] = inimoms[1]
                    inimoms[1] = inimoms[2]
                    inimoms[2] = moment
                if moment == 0:
                    moment = np.mean(inimoms[np.nonzero(inimoms)])
            
            if phs == 'P':
                outpfile.write('%.8e %s\n' % (moment,sta))
            elif phs == 'S' and dstime[sta] != -1:
                outsfile.write('%.8e %s\n' % (moment,sta))

    outpfile.close()
    outsfile.close()