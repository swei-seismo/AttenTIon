import os
import momcalc
import subprocess
import numpy as np
from mpi4py import MPI
from pathlib import Path


def load_catalog(catalog_file):
    catalog_dict = {}
    with open(catalog_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            orid = parts[-1]
            catalog_dict[orid] = {
                'evlo': float(parts[1]),
                'evla': float(parts[2]),
                'evdp': float(parts[3])
            }

    return catalog_dict


def get_station_data(sacdir, orid):
    cmd = f"saclst KSTNM t1 gcarc f {sacdir}/{orid}/*Z.sac"
    output = subprocess.check_output(cmd, shell=True, text=True).strip()
    dbev = {}
    dstime = {}
    for line in output.split("\n"):
        parts = line.strip().split()
        sta = parts[1]
        stime = parts[2]
        gcarc = parts[3]
        dbev[sta] = float(gcarc)
        dstime[sta] = float(stime)
    
    return dbev, dstime


def calculate_gs_for_event(orid, catalog_dict, sacdir, workdir):
    depth = catalog_dict[orid]['evdp']
    depth = max(depth, 0.00)
    dbev, dstime = get_station_data(sacdir, orid)
    sorted_stations = sorted(dbev.items(), key=lambda item: item[1])
    num = len(sorted_stations)
    inimomp = np.array([1.0e-12, 5.0e-13, 1.0e-13])
    inimoms = np.array([1.0e-12, 5.0e-13, 1.0e-13])
    p_data = []
    s_data = []
    for i in range(num):
        sta = sorted_stations[i][0]
        delta = sorted_stations[i][1]
        moment_p = momcalc.correct(delta, depth, 'P')
        if moment_p != 0:
            inimomp[0] = inimomp[1]
            inimomp[1] = inimomp[2]
            inimomp[2] = moment_p
        else:
            moment_p = np.mean(inimomp[np.nonzero(inimomp)])
        
        p_data.append(f'{moment_p:.8e} {sta}\n')
        
        if dstime[sta] != -12345:
            moment_s = momcalc.correct(delta, depth, 'S')
            if moment_s != 0:
                inimoms[0] = inimoms[1]
                inimoms[1] = inimoms[2]
                inimoms[2] = moment_s
            else:
                moment_s = np.mean(inimoms[np.nonzero(inimoms)])
            
            if dstime[sta] != -1:
                s_data.append(f'{moment_s:.8e} {sta}\n')
    
    outpfln = Path(workdir) / f'pgsfile_{orid}.txt'
    outsfln = Path(workdir) / f'sgsfile_{orid}.txt'
    
    with open(outpfln, 'w') as f:
        f.writelines(p_data)
    
    with open(outsfln, 'w') as f:
        f.writelines(s_data)


print("Calculating constant")
main_path = Path.cwd()
workdir = main_path / "GS"
namedir = main_path / "extra_Eventname"
catalog_file = main_path / "AACSE_catalog.dat"
sacdir = main_path / "correctedSeismograms"
workdir.mkdir(parents=True, exist_ok=True)
print("Loading catalog...")
catalog_dict = load_catalog(catalog_file)
print(f"Loaded {len(catalog_dict)} events from catalog")
with open(namedir, 'r') as f:
    oridlst = [line.strip() for line in f if line.strip()]

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

events_to_process = [oridlst[i] for i in range(len(oridlst)) if i % size == rank]

print(f"Rank {rank}: Processing {len(events_to_process)} events out of {len(oridlst)} total")
for idx, orid in enumerate(events_to_process):
    if idx % 10 == 0:
        print(f"Rank {rank}: Processing event {idx+1}/{len(events_to_process)}: {orid}")
    
    calculate_gs_for_event(orid, catalog_dict, sacdir, workdir)

print(f"Rank {rank}: Completed processing {len(events_to_process)} events")