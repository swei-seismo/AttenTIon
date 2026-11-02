import os
import glob
import datetime
import subprocess
import numpy as np
from mpi4py import MPI
from pathlib import Path


def load_catalog(catalog_file):
    catalog_dict = {}
    with open(catalog_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            event_id = parts[-1]
            catalog_dict[event_id] = {
                'evlo': float(parts[1]),
                'evla': float(parts[2]),
                'evdp': float(parts[3])
            }
    
    return catalog_dict


def load_picks(picks_file):
    """Load picks into nested dictionary for fast lookup"""
    picks_dict = {}
    with open(picks_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            orid = parts[0]
            sta = parts[2]
            chan = parts[3]
            phase = parts[6]
            
            if orid not in picks_dict:
                picks_dict[orid] = {}
            if sta not in picks_dict[orid]:
                picks_dict[orid][sta] = []
            
            pick_data = {
                'origintime': parts[4],
                'picktime': parts[5],
                'chan': chan,
                'phase': phase
            }
            picks_dict[orid][sta].append(pick_data)
    
    return picks_dict


def check_single_channel():
    cmd = "ls *Z.sac | awk -F '.' '{print $1, $2}'"
    output = subprocess.check_output(cmd, shell=True, text=True).strip()
    lines = output.split("\n")
    sac_commands = ""
    for line in lines:
        line = line.strip()
        parts = line.split()
        net = parts[0]
        sta = parts[1]
        lst = glob.glob(f"{net}.{sta}.*.sac")
        if len(lst) == 3:
            continue

        if len(lst) == 2:
            print(net, sta)
        
        for sacfile in lst:
            chan = sacfile.split(".")[2]
            new_chan1 = chan[:-1]+"E"
            new_sacfile1 = f"{net}.{sta}.{new_chan1}.sac"
            sac_commands += f"read {sacfile}\n"
            sac_commands += f"ch KCMPNM {new_chan1}\n"
            sac_commands += f"w {new_sacfile1}\n"

            new_chan2 = chan[:-1]+"N"
            new_sacfile2 = f"{net}.{sta}.{new_chan2}.sac"
            sac_commands += f"read {sacfile}\n"
            sac_commands += f"ch KCMPNM {new_chan2}\n"
            sac_commands += f"w {new_sacfile2}\n"
    
    if sac_commands:
        sac_commands += "q \n"
        subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(sac_commands.encode())


def add_sac_headers(sac_z_files, event_id, catalog_dict, picks_dict):

    event_data = catalog_dict[event_id]
    evlo = event_data['evlo']
    evla = event_data['evla']
    evdp = event_data['evdp']
    picks_for_event = picks_dict.get(event_id, {})
    sac_commands = ""
    for sacfile in sac_z_files:
        net = sacfile.split(".")[0]
        sta = sacfile.split(".")[1]
        sac_commands += f"readhdr {net}.{sta}.*.sac\n"
        sac_commands += f"ch evla {evla} evlo {evlo} evdp {evdp}\n"
        sac_commands += "ch lcalda True\n"
        if sta in picks_for_event:
            picks_list = picks_for_event[sta]
            origin_time_str = None
            p_pick = None
            s_pick = None

            for pick in picks_list:
                if origin_time_str is None:
                    origin_time_str = pick['origintime']
                
                if pick['phase'] == 'P' and 'Z' in pick['chan']:
                    p_pick = pick
                
                if pick['phase'] == 'S' and s_pick is None:
                    s_pick = pick
            
            if origin_time_str:
                o = datetime.datetime.strptime(origin_time_str, "%Y-%m-%dT%H:%M:%S.%fZ")
                o_jday = o.strftime("%j")
                o_msec = int(o.microsecond / 1000 + 0.5)
                sac_commands += f"ch o gmt {o.year} {o_jday} {o.hour} {o.minute} {o.second} {o_msec}\n"
            
            # P arrival (T0)
            if p_pick:
                t0 = datetime.datetime.strptime(p_pick['picktime'], "%Y-%m-%dT%H:%M:%S.%fZ")
                t0_jday = t0.strftime("%j")
                t0_msec = int(t0.microsecond / 1000 + 0.5)
                sac_commands += f"ch t0 gmt {t0.year} {t0_jday} {t0.hour} {t0.minute} {t0.second} {t0_msec}\n"
            
            # S arrival (T1)
            if s_pick:
                t1 = datetime.datetime.strptime(s_pick['picktime'], "%Y-%m-%dT%H:%M:%S.%fZ")
                t1_jday = t1.strftime("%j")
                t1_msec = int(t1.microsecond / 1000 + 0.5)
                sac_commands += f"ch t1 gmt {t1.year} {t1_jday} {t1.hour} {t1.minute} {t1.second} {t1_msec}\n"
            
        sac_commands += "wh \n"
    
    if sac_commands:
        sac_commands += "q \n"
        subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(sac_commands.encode())


def check_and_remove_no_P_arrival():
    cmd = "saclst t0 f *.sac"
    output = subprocess.check_output(cmd, shell=True, text=True).strip()
    for line in output.split('\n'):
        parts = line.split()
        filename = parts[0]
        t0_value = float(parts[1])
        if t0_value == -12345:
            os.remove(filename)


def sactotxt(sacfiles):
    sac_commands = ""
    for sacfile in sacfiles:
        txt_filename = sacfile.replace("sac", "txt")
        sac_commands += f"r {sacfile}\n"
        sac_commands += "ch allt (0 - &1,T0&) iztype IT0\n"
        sac_commands += f"w alpha {txt_filename}\n"
    
    if sac_commands:
        sac_commands += "q\n"
        subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(sac_commands.encode())


def process_event(dirname, seismograms, catalog_dict, picks_dict):
    event_path = Path(seismograms)/dirname
    os.chdir(event_path)
    
    check_single_channel()

    sac_z_files = glob.glob("*Z.sac")
    if sac_z_files:
        add_sac_headers(sac_z_files, dirname, catalog_dict, picks_dict)
    
    check_and_remove_no_P_arrival()

    all_sacfiles = glob.glob("*.sac")
    if all_sacfiles:
        sactotxt(all_sacfiles)


os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
main_path = os.getcwd()
catalog_file = Path(main_path)/"AACSE_catalog.dat"
picks_file = Path(main_path)/"new_picks.dat"
seismograms = Path(main_path)/"correctedSeismograms"
catalog_dict = load_catalog(catalog_file)
picks_dict = load_picks(picks_file)
print(f"Loaded {len(catalog_dict)} events and picks for {sum(len(v) for v in picks_dict.values())} stations")
with open('extra_Eventname', 'r') as f:
    oridlst = [line.strip() for line in f if line.strip()]

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
events_to_process = [oridlst[i] for i in range(len(oridlst)) if i % size == rank]
print(f"Rank {rank}: Processing {len(events_to_process)} events out of {len(oridlst)} total")
for idx, dirname in enumerate(events_to_process):
    if idx % 10 == 0:
        print(f"Rank {rank}: Processing event {idx+1}/{len(events_to_process)}: {dirname}")
    
    try:
        process_event(dirname, seismograms, catalog_dict, picks_dict)
    except Exception as e:
        print(f"Rank {rank}: Error processing {dirname}: {e}")
        continue

print(f"Rank {rank}: Completed processing {len(events_to_process)} events")