# Program to process sacfiles

import os
import glob
import datetime
import subprocess
import numpy as np
from mpi4py import MPI


def check_single_channel():
    ## There is only 1 channel(EHZ) in some stations, 
    ## Copy vertical channel and make them into two horizontal channels
    cmd = "ls *Z.sac | awk -F '.' '{print $1, $2}'"
    output = os.popen(cmd).read().strip().split('\n')
    for line in output:
        line = line.strip()
        net = line.split()[0]
        sta = line.split()[1]
        lst = glob.glob("%s.%s.*.sac" %(net, sta))
        if len(lst) == 3:
            continue
        else:
            for sacfile in lst:
                chan = sacfile.split(".")[2]
                new_chan1 = chan[:-1]+"E"
                new_sacfile1 = net+"."+sta+"."+new_chan1+".sac"
                s = ""
                s += "read %s \n" %(sacfile)
                s += "ch KCMPNM %s \n" %(new_chan1)
                s += "w %s \n" %(new_sacfile1)
                s += "q \n"
                subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(s.encode())

                new_chan2 = chan[:-1]+"N"
                new_sacfile2 = net+"."+sta+"."+new_chan2+".sac"
                s = ""
                s += "read %s \n" %(sacfile)
                s += "ch KCMPNM %s \n" %(new_chan2)
                s += "w %s \n" %(new_sacfile2)
                s += "q \n"
                subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(s.encode())


def add_sac_header(sacfile, event_id, catalog, picks_file):
    net = sacfile.split(".")[0]
    sta = sacfile.split(".")[1]

    # evlo, evla, evdp
    with open(catalog, "r") as f:
        for line in f:
            if event_id in line:
                parts = line.strip().split()
                evlo = float(parts[1])
                evla = float(parts[2])
                evdp = float(parts[3])
                break
    

    # Origintime (o), P(T0) and S(T1) wave arrival
    o = 0
    t0 = 0
    t1 = 0
    with open(picks_file, "r") as f:
        for line in f:
            ## You need to check if there is two stations like "AUL" and "AULG"
            ## Some errors could happen due to "sta in line"
            if event_id in line and sta in line:
                parts = line.strip().split()
                phase = parts[6]
                origintime = parts[4]
                picktime = parts[5]
                chan = parts[3]
                o = datetime.datetime.strptime(origintime, "%Y-%m-%dT%H:%M:%S.%fZ")
                o_jday = o.strftime("%j")
                o_msec = int(o.microsecond/1000+0.5)
                if phase == "P" and "Z" in chan:
                    t0 = datetime.datetime.strptime(picktime, "%Y-%m-%dT%H:%M:%S.%fZ")
                    t0_jday = t0.strftime("%j")
                    t0_msec = int(t0.microsecond/1000+0.5)
                elif phase == "S":
                    t1 = datetime.datetime.strptime(picktime, "%Y-%m-%dT%H:%M:%S.%fZ")
                    t1_jday = t1.strftime("%j")
                    t1_msec = int(t1.microsecond/1000+0.5)
    
    s = ""
    s += "readhdr %s.%s.*.sac \n" %(net,sta)
    s += "ch evla %s \n" %(evla)
    s += "ch evlo %s \n" %(evlo)
    s += "ch evdp %s \n" %(evdp)
    s += "ch lcalda True \n" # calculate distance and azimuth
    if type(o) is datetime.datetime:
        s += "ch o gmt {} {} {} {} {} {} \n".format(o.year,o_jday,o.hour,
                                             o.minute,o.second,o_msec)
    if type(t0) is datetime.datetime:
        s += "ch t0 gmt {} {} {} {} {} {} \n".format(t0.year,t0_jday,t0.hour,
                                              t0.minute,t0.second,t0_msec)
    if type(t1) is datetime.datetime:
        s += "ch t1 gmt {} {} {} {} {} {} \n".format(t1.year,t1_jday,t1.hour,
                                              t1.minute,t1.second,t1_msec)
    s += "wh \n"
    s += "q \n"
    subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.DEVNULL).communicate(s.encode())


def check_T0(sacfile):
    ## Some stations only have S arrival without P arrival
    ## I don't use them due to following Qp/Qs ratio measurements
    command = "saclst t0 f %s" %(sacfile)
    output = os.popen(command).read().strip()
    T0 = float(output.split()[1])
    if T0 == -12345:
        os.remove("./%s" %(sacfile))


def sactotxt(sacfile):
    txt_filename = sacfile.replace("sac","txt")
    s = ""
    s += "r %s \n" %(sacfile)
    s += "ch allt (0 - &1,T0&) iztype IT0 \n" # set T0 as zero time
    s += "w alpha %s \n" %(txt_filename)
    s += "q \n"
    subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.DEVNULL).communicate(s.encode())


os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
main_path = os.getcwd()
catalog = main_path+"/AACSE_catalog.dat"
picks_file = main_path+"/picks.dat"
seismograms = main_path+"/correctedSeismograms"
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

for idir in np.arange(bgn, bgn+NUM, 1):
    dirname = oridlst[idir]
    print(f"Rank {rank}: {dirname}")
    os.chdir("%s/%s" %(seismograms, dirname))

    check_single_channel()

    sac_list = glob.glob("*Z.sac")
    for sac in sac_list:
        add_sac_header(sac,dirname,catalog,picks_file)
    
    sacfls = glob.glob("*.sac")
    for sacfile in sacfls:
        check_T0(sacfile)
    
    sacfls = glob.glob("*.sac")    
    for sacfile in sacfls:
        sactotxt(sacfile)