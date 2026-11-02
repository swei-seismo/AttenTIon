import os
import glob
import subprocess
import numpy as np
from mpi4py import MPI
from pathlib import Path

def load_pz_metadata(PZfile_path):
    pz_metadata = {}
    pz_files = list(Path(PZfile_path).glob("SAC_PZ_*"))
    for pz_file in pz_files:
        parts = pz_file.stem.split("_")
        net = parts[2]
        sta = parts[3]
        key = f"{net}.{sta}"
        
        with open(pz_file, "r") as f:
            content = f.read()
            for line in content.split("\n"):
                if 'LATITUDE' in line and ':' in line:
                    stla = float(line.split(':')[-1].strip())
                elif 'LONGITUDE' in line and ':' in line:
                    stlo = float(line.split(':')[-1].strip())
            
            pz_metadata[key] = {
                'file': str(pz_file),
                'stla': stla,
                'stlo': stlo
            }
    print(f"Loaded metadata for {len(pz_metadata)} PZ files")

    return pz_metadata


def remove_instrument_response_batch(sac_files, PZfile_path, pz_metadata, output_dir):
    files_str = " ".join(f.name for f in sac_files)
    print(files_str)
    cmd = f"saclst delta f {files_str}"
    output = subprocess.check_output(cmd, shell=True, text=True).strip()
    delta_dict = {}
    for line in output.split("\n"):
        parts = line.split()
        delta_dict[parts[0]] = float(parts[1])
    
    sac_commands = ""
    files_processed = 0
    for sac_path_obj in sac_files:
        sacfile = sac_path_obj.name
        delta = delta_dict.get(sacfile)
        sampling_rate = 1.0/delta

        basename = Path(sacfile).name.replace("SAC", "sac")
        new_sacname = Path(output_dir) / basename
        parts = sacfile.split(".")
        net = parts[0]
        sta = parts[1]
        key = f"{net}.{sta}"
        if key not in pz_metadata:
            print(f"  Skipping {sacfile}: PZ file not found for {key}")
            continue
        pz_info = pz_metadata[key]
        pz_file = pz_info['file']
        stla = pz_info['stla']
        stlo = pz_info['stlo']
        
        if 38 <= sampling_rate <= 42:  # ~40 Hz
            freq_params = "0.01 0.05 18 19"
        elif 48 <= sampling_rate <= 52:  # ~50 Hz
            freq_params = "0.01 0.05 22 24"
        elif 98 <= sampling_rate <= 102:  # ~100 Hz
            freq_params = "0.01 0.05 45 48"
        else:
            print(f"  Skipping {sacfile}: unsupported sampling rate {sampling_rate:.1f} Hz")
            continue
        
        sac_commands += f"r {sacfile}\n"
        sac_commands += "rmean; rtr; taper\n"
        sac_commands += f"trans from polezero subtype {pz_file} to vel freq {freq_params}\n"
        sac_commands += "mul 1.0e9\n"
        sac_commands += f"ch stla {stla}\n"
        sac_commands += f"ch stlo {stlo}\n"
        sac_commands += f"w {new_sacname}\n"
        
        files_processed += 1
        
    if sac_commands:
        sac_commands += "q\n"
        subprocess.Popen(['sac'], 
                        stdin=subprocess.PIPE,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL).communicate(sac_commands.encode())
    
    print(f"  Processed {files_processed} files")


print("Removing instrument response")
os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
main_path = Path.cwd()
sac_path = main_path / "processedSeismograms"
output_base = main_path / "correctedSeismograms"
PZfile_path = Path("/mnt/research/seismo_wei/ZZhang/tstar_inversion/tstar_data/SACPZ")
print("Loading PZ file metadata...")
pz_metadata = load_pz_metadata(PZfile_path)
with open('extra_Eventname', 'r') as f:
    oridlst = [line.strip() for line in f if line.strip()]

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

events_to_process = [oridlst[i] for i in range(len(oridlst)) if i % size == rank]
print(f"Rank {rank}: Processing {len(events_to_process)} events out of {len(oridlst)} total")

for idx, dirname in enumerate(events_to_process):
    input_dir = sac_path / dirname
    output_dir = output_base / dirname
    if idx % 10 == 0:
        print(f"Rank {rank}: Processing event {idx+1}/{len(events_to_process)}: {dirname}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(input_dir)
    sac_list = list(input_dir.glob("*.SAC"))
    remove_instrument_response_batch(sac_list, PZfile_path, pz_metadata, output_dir)

print(f"Rank {rank}: Completed processing {len(events_to_process)} events")
        