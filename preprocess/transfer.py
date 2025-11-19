# Program to remove instrument response

import os
import glob
import subprocess
import numpy as np
from mpi4py import MPI


def remove_instrument_response(sacfile, PZfile_path, output_dir):
    cmd = "saclst delta f %s | awk -F ' ' '{print $2}'" %(sacfile)
    sampling_rate = 1/float(os.popen(cmd).read().strip())

    basename = os.path.basename(sacfile).replace("SAC", "sac")
    new_sacname = os.path.join(output_dir, basename)

    net = sacfile.split(".")[0]  
    sta = sacfile.split(".")[1]
    pz = glob.glob("%s/SAC_PZ_%s_%s" %(PZfile_path, net, sta))
    if len(pz) != 1:
        print("PZ file error for %s" %(sacfile))
        return

    cmd2 = "cat %s | awk '/LATITUDE/ {print}'" %(pz[0])
    output2 = os.popen(cmd2).read().strip()
    stla = float(output2.split(":")[-1].strip())

    cmd3 = "cat %s | awk '/LONGITUDE/ {print}'" %(pz[0])
    output3 = os.popen(cmd3).read().strip()
    stlo = float(output3.split(":")[-1].strip())

    s = ""
    s += "r %s \n" %(sacfile)
    s += "rmean; rtr; taper \n"
    if sampling_rate == 40:
        s += "trans from polezero subtype %s to vel freq 0.01 0.05 18 19 \n" %(pz[0])
    elif sampling_rate == 50:
        s += "trans from polezero subtype %s to vel freq 0.01 0.05 22 24 \n" %(pz[0])
    elif sampling_rate == 100:
        s += "trans from polezero subtype %s to vel freq 0.01 0.05 45 48 \n" %(pz[0])
    s += "mul 1.0e9 \n"
    s += "ch stla %s \n" %(stla)
    s += "ch stlo %s \n" %(stlo)
    s += "w %s \n" %(new_sacname)
    s += "q \n"
    subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.DEVNULL).communicate(s.encode())


os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
main_path = os.getcwd()
sac_path = main_path+"/processedSeismograms"
output_base = os.path.join(main_path, "correctedSeismograms")
PZfile_path = "/mnt/research/seismo_wei/ZZhang/tstar_inversion/tstar_data/SACPZ"
oridlst = []
forid = open('Eventname', 'r')
for name in forid.readlines():
    name = name.strip('\n')
    if name is not None:
        oridlst.append(name)
forid.close()

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
num = int(len(oridlst)/size)

bgn = rank*num
num_last = len(oridlst)-(size-1)*num
NUM = (rank<size-1 and num or num_last)

for idir in np.arange(bgn,bgn+NUM,1):
    dirname = oridlst[idir]
    print(f"Rank {rank}: {dirname}")

    input_dir = os.path.join(sac_path, dirname)
    output_dir = os.path.join(output_base, dirname)
    os.makedirs(output_dir, exist_ok=True)

    os.chdir(input_dir)
    sac_list = glob.glob("*.SAC")
    for sacfile in sac_list:
        remove_instrument_response(sacfile, PZfile_path, output_dir)