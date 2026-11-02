import os
import glob
import numpy as np
from mpi4py import MPI
from pathlib import Path
from obspy import read, Stream


def mseedtosac(miniseed_files, sac_path):
    Path(sac_path).mkdir(parents=True, exist_ok=True)
    for miniseed in miniseed_files:
        st = read(miniseed, headonly=False)
        tr = st[0]
        net = tr.stats.network
        sta = tr.stats.station
        cha = tr.stats.channel
        sacfile = f"{sac_path}/{net}.{sta}.{cha}.SAC"
        single_trace = Stream(traces=[tr])
        single_trace.write(sacfile, format="SAC")


main_path = Path.cwd()
miniseed_path = Path("/mnt/research/seismo_wei/ZZhang/tstar_inversion/tstar_data/waveforms")
with open('extra_Eventname', 'r') as f:
    oridlst = [line.strip() for line in f if line.strip()]
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    events_to_process = [oridlst[i] for i in range(len(oridlst)) if i % size == rank]
    
    print(f"Rank {rank}: Processing {len(events_to_process)} events out of {len(oridlst)} total")
    
    for idx, dirname in enumerate(events_to_process):
        event_path = miniseed_path / dirname
        
        if idx % 10 == 0:
            print(f"Rank {rank}: Processing event {idx+1}/{len(events_to_process)}: {dirname}")
        
        os.chdir(event_path)
        miniseed_list = list(event_path.glob("*.mseed"))
        sac_path = main_path / "processedSeismograms" / dirname
        mseedtosac(miniseed_list, sac_path)

print(f"Rank {rank}: Completed processing {len(events_to_process)} events")
