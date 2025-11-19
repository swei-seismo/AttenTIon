import os
import glob
import numpy as np
from obspy import read
from mpi4py import MPI


def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def mseedtosac(miniseed, sac_path):
    st = read(miniseed)
    tr = st[0]
    net = tr.stats.network
    sta = tr.stats.station
    cha = tr.stats.channel
    sacfile = "%s/%s.%s.%s.SAC" %(sac_path, net, sta, cha)
    st.write(sacfile, format="SAC")


main_path = os.getcwd()
miniseed_path = "/mnt/research/seismo_wei/ZZhang/tstar_inversion/tstar_data/waveforms"

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
    os.chdir("%s/%s" %(miniseed_path,dirname))
    print(f"Rank {rank}: {dirname}")
    sac_path = main_path+"/processedSeismograms/"+dirname
    mkdir(sac_path)
    miniseed_list = glob.glob("*.mseed")
    for miniseed in miniseed_list:
        mseedtosac(miniseed, sac_path)
