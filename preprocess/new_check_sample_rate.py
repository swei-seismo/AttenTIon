import os
import subprocess
import numpy as np
from mpi4py import MPI
from pathlib import Path


def check_and_fix_sample_rate():
    cmd = "saclst delta f *.SAC"
    output = subprocess.check_output(cmd, shell=True, text=True).strip()
    lines = output.splitlines()
    possible_deltas = {0.01, 0.02, 0.025}

    files_to_fix = []
    for line in lines:
        parts = line.split()
        sacname = parts[0]
        delta = float(parts[1])
        if delta not in possible_deltas:
            print(f"  Fixing {sacname}: delta={delta} -> 0.025")
            files_to_fix.append(sacname)
    
    if files_to_fix:
        sac_commands = ""
        for sacfile in files_to_fix:
            sac_commands += f"r {sacfile}\n"
            sac_commands += "ch delta 0.025\n"
            sac_commands += f"w {sacfile}\n"
        
        sac_commands += "q \n"
        subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(sac_commands.encode()) 


print("Checking and fixing sample rates")
os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
main_path = Path.cwd()
sac_path = main_path / "processedSeismograms"
with open('extra_Eventname', 'r') as f:
    oridlst = [line.strip() for line in f if line.strip()]

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

events_to_process = [oridlst[i] for i in range(len(oridlst)) if i % size == rank]

print(f"Rank {rank}: Processing {len(events_to_process)} events out of {len(oridlst)} total")

for idx, dirname in enumerate(events_to_process):
    event_path = sac_path / dirname
    
    if idx % 10 == 0:
        print(f"Rank {rank}: Processing event {idx+1}/{len(events_to_process)}: {dirname}")
    
    os.chdir(event_path)
    check_and_fix_sample_rate()

print(f"Rank {rank}: Completed processing {len(events_to_process)} events")