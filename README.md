# t* inversion program by MSU seismo

## History (some early history might be missing):
2002, Josh Stachnik wrote a MATLAB package to invert for t* (Stachnick et al., 2004, JGR)

2007, Sarah Pozgay made significant updates on the MATLAB package (Pozgay et al., 2009, G-cubed)

2013-2015, Shawn Wei rewrote the MATLAB package in Python and made significant changes (Wei & Wiens, 2018, EPSL; Wei & Wiens, 2020, JGR)

2019, Steve Hicks made changes (Hicks et al., in prep.)

2020-2022, Yurong Zhang and Zhuoran Zhang ....


# t* inversion

### Step 1. Pre_process

Usually we will have two formats of data used in this inversion, which end with ".sac" and ".txt". The sac files should contain the sac headers as below:

1. KNETWK (network name)
2. KSTNM  (station name)
3. KCMPNM (channel name)
4. stla   (latitude of station)
5. stlo   (longitude of station)
6. evla   (latitude of event)
7. evlo   (longitude of event)
8. evdp   (depth of event)
9. dist   (distance in kilometers)
10. gcarc (distance in degree)
11. baz   (back azimuth)
12. delta (sampling rate)
13. o     (origin time)
14. T0    (P wave arrival time)
15. T1    (S wave arrival time)

The processed seismograms are stored in `./data/processedSeismograms`

### Step 2. Run mkgsfl.py (Subroutines: momcalc.py)

Calculate geometrical spreading and free surface effects for P and S waves. The files are named as `pgsfile_eventid` and `sgsfile_eventid`

The calculated files are stored in `./data/GS`

### Step 3. Run main_tstar_inversion.py (Subroutines: tstar_parameters.py, tstar_load.py, tstarsub.py, tstar_inversion_function.py)

OUTPUT:
(eventfocal027.log)  ==> corner frequency and magnitude for each event.

(bestfc.lst)         ==> corner frequency for P wave, magnitude and ln(M0) (Seismic moment) for each event.

(result)             ==> spectrum, t* and misfit for each event-station pair.

(plotfall)           ==> corner frequency versus L2 norm and log10(moment) for each event.

(plottstar-fc)       ==> corner frequency versus t* perturbation for each event-station pair.

(plotspec)           ==> spectrum in frequency domain.

(plotseis)           ==> seismogram in time domain.

(plotsnr)            ==> signal-to-noise ratio of spectrum.

