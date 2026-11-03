# AttenTIon t* inversion workflow

This repository contains the preprocessing and fixed-corner-frequency inversion used to estimate P- and S-wave attenuation factor (`t*`) of the event-station pairs in the Alaska Peninsula of the Alaska-Aluetian Subduction Zone. 

The preprocessing scripts (1) download the waveform (`new_download.py`) (2) convert miniseed files to SAC files (`new_mseed2sac.py`), (3) Remove the instrument response (`new_transfer.py`), (4) write timestamp into header of SAC files and convert SAC files into txt format (`new_process_sac.py`) .


The scripts currently assume the research-space layout under `/mnt/research/seismo_wei/ZZhang/tstar_inversion`. In particular, most waveform and input paths are set in `tstar_parameters.py`.

## Software and working directories

- Python packages: `numpy`, `scipy`, `matplotlib`, `mpi4py`, `obspy` `pandas` (imported by `tstar_load.py`), and mtspec (https://krischer.github.io/mtspec/)
- External commands: SAC (`sac` and `saclst`) and an MPI launcher such as `srun`.
- The inversion imports the bundled `lib/mtspec*.so` through `multitaper.py`. This compiled library and its dependencies must be compatible with the Python environment and machine used to run the code. It usually works well when using different versions of Python.
- Run `preprocess/new*.py` with the data directory as the **current working directory**. Run `main.py` with an inversion directory containing `Eventname` as the **current working directory**. Inversion outputs are written beneath that directory's `workdir/`.

## Inputs 
#### The input files demonstrate the expected format and contain only a subset of the full dataset to reduce the repository size.
| Input | Expected content and use |
| --- | --- |
| `Eventname` | One event ID per line; selects events for the `new*` preprocessing scripts and processed by `main.py` |
| `download_data_list.dat` | Event ID, station, and UTC origin time per line; `new_download.py` uses it to select stations and request waveforms from 60 s before to 300 s after the origin. |
| `AACSE_catalog.dat` | Origin time, longitude, latitude, depth (km), magnitude, event ID; used for event headers, geometry, and inversion event information. |
| `new_picks.dat` | Event ID, network, station, channel, origin time, pick time, phase (`P` or `S`), travel time; used to set SAC `o`, `t0` (P), and `t1` (S). |
| `SACPZ/SAC_PZ_*` | Station pole-zero responses and station coordinates for response removal. |
| `Alaska_events_source_spectra.dat` | Source information keyed by event ID; column 6 supplies the independently constrained corner frequency.

The default preprocessing data directory is `/mnt/research/seismo_wei/ZZhang/tstar_inversion/tstar_data`. `new_mseed2sac.py` and `new_transfer.py` also contain absolute paths to its `waveforms/` and `SACPZ/` directories, respectively. The inversion expects `correctedSeismograms/`, `GS/`, `AACSE_catalog.dat`, and `Alaska_events_source_spectra.dat` under the same data directory; edit `tstar_parameters.py` to use another location.


## Preprocessing sequence

Run the scripts in this order from the data directory. The `new_download.py` step is serial; the later scripts distribute event IDs across MPI ranks.

| Step | Script | Main output |
| --- | --- | --- |
| 1 | `preprocess/new_download.py` | `waveforms/<event_id>/*.mseed` and station metadata in `stations/`. Its configured request region is 51–60°N, 163–147°W. |
| 2 | `preprocess/new_mseed2sac.py` | `processedSeismograms/<event_id>/<network>.<station>.<channel>.SAC`. |
| 3 | `preprocess/new_check_sample_rate.py` | Updates the SAC `delta` header in `processedSeismograms/` to 0.025 s when it is not 0.025 s. I have some traces with sample rate of 0.025001s, so I need to correct this problem. |
| 4 | `preprocess/new_transfer.py` | Removes instrument response using SAC pole-zero files, converts to velocity, multiplies by `1e9`, and writes `correctedSeismograms/<event_id>/*.sac`. Traces lacking a matching pole-zero file or a supported sample rate are skipped. |
| 5 | `preprocess/new_process_sac.py` | Adds event coordinates and origin/P/S times to SAC headers, removes traces without a P `t0` header, and writes corresponding SAC alpha `*.txt` files used by the inversion. |
| 6 | `preprocess/new_mkgsfl.py` | Writes `GS/pgsfile_<event_id>.txt` and `GS/sgsfile_<event_id>.txt`, each containing a correction factor and station name. `momcalc.py` computes these factors from the iasp91 model, source depth, and event–station distance. |

For example, inside a suitable Slurm allocation, replace `/path/to/AttenTIon` with the clone location:

```bash
cd /mnt/research/seismo_wei/ZZhang/tstar_inversion/tstar_data
python /path/to/AttenTIon/preprocess/new_download.py
srun -n 20 python /path/to/AttenTIon/preprocess/new_mseed2sac.py
srun -n 20 python /path/to/AttenTIon/preprocess/new_check_sample_rate.py
srun -n 20 python /path/to/AttenTIon/preprocess/new_transfer.py
srun -n 20 python /path/to/AttenTIon/preprocess/new_process_sac.py
srun -n 20 python /path/to/AttenTIon/preprocess/new_mkgsfl.py
```

`new_process_sac.py` may create E/N files by copying vertical channel data, in this case, the S-wave t* measurement comes from the S waveform on the vertical channel. This script will also remove of traces without P picks changes files in `correctedSeismograms/`.


## Inversion strategy

The active settings are in `tstar_parameters.py`: `source_para = 1`, fixed `alpha = 0.27`, `fcps = 1.0`, and `add_site = 0`. Despite the older parameter comment about a corner-frequency grid search, this workflow reads the event's corner frequency from column 6 of `Alaska_events_source_spectra.dat` and uses that value. With `fcps = 1.0`, P and S use the same corner frequency. No separate site-spectrum correction is applied.

1. `main.py` loads the event, SAC arrivals, three-component text traces, and P/S geometric-spreading and free-surface factors. P uses the vertical component; S is selected from the horizontal components. A P pick is required, but an S pick is optional.
2. `tstarsub.py` takes 3 s P and 6 s S signal windows with preceding noise windows. It calculates sine-multitaper spectra, converts velocity spectra to displacement amplitude, and selects a contiguous frequency band using spectral SNR, bandwidth, and negative-slope/correlation tests. The configured SNR thresholds are 3 for P and 2 for S.
3. For a fixed corner frequency, the modeled amplitude is proportional to `C × M0 / (1 + (f/fc)^2) × exp[-π f^(1−alpha) t*]`. Here `C` includes the geometric-spreading/free-surface factor, `M0` is event moment, and `t*` is estimated for each event–station path. The inversion solves the logarithmic form, in which `ln(M0)` and path `t*` enter linearly.
4. A preliminary bounded least-squares fit includes an event–station constant log-amplitude offset for each path. Stations with a spectral-fit score below 0.75 are removed. At least six accepted P stations are needed to proceed; S is attempted after P succeeds and needs at least three accepted S stations.
5. The constant offsets are re-estimated and folded into updated amplitude corrections. The final fit holds those corrections fixed and uses non-negative least squares for `ln(M0)` and path `t*`, separately for P and S. The script also computes formal errors, residual diagnostics, moment magnitude, and the reported `Qp/Qs` ratio from the P and S attenuation estimates.

## Inversion outputs

Run `main.py` from a directory with `Eventname`. The code creates `workdir/` there. A typical Slurm command inside an allocation is:

```bash
cd /mnt/research/seismo_wei/ZZhang/tstar_inversion/fixed_fc_with_constant_correction
srun -n 20 python /path/to/AttenTIon/main.py
```

With the default `alpha = 0.27`, output filenames end in `027`:

| File | Contents |
| --- | --- |
| `workdir/result/<event_id>_pstar027.dat` | One row per accepted P station: event ID, station, latitude, longitude, depth, P `t*`, formal error, misfit, and `1000 × t*_P / P_travel_time`. |
| `workdir/result/<event_id>_sstar027.dat` | One row per accepted S station: event ID, station, location/depth, S `t*`/error/misfit/attenuation, P counterparts, and the reported `Qp/Qs` ratio. |
| `workdir/result/<event_id>_perr027.dat`, `*_serr027.dat` | Residual norm, number of spectral data points, and two normalized residual measures for the final P or S fit. |
| `workdir/fitting.lst` | Event ID, station, phase, and preliminary spectral-fit score. |
| `workdir/error.lst` | Event ID, station, phase, original log correction, and estimated constant offset; the active write is for P. |
| `workdir/event.log` | Per-event progress and P/S moment-magnitude information. |

Results are overwritten per event, while `event.log`, `fitting.lst`, and `error.lst` are opened in append mode. Plotting is disabled by default; the plotting flags in `tstar_parameters.py` can enable figures under `workdir/specfig/`.
