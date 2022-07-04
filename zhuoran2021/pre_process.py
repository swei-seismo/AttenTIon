import os
import glob
import subprocess
import datetime
from obspy import read
from obspy import read_inventory


def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def mseedtosac(miniseed,xmlfile_path,sac_path):
    # Convert miniseed file to SAC file
    # Remove instrument response
    # Frequency band of interest: 100s - 50Hz
    original_st = read(miniseed)
    st = original_st.copy()
    tr = st[0]
    net = tr.stats.network
    sta = tr.stats.station
    loc = tr.stats.location
    cha = tr.stats.channel
    xmlfile = xmlfile_path+"/%s.%s.xml" %(net,sta)
    inv = read_inventory(xmlfile).select(location=loc, channel=cha)
    try:
        print(net,sta,xmlfile)
        tr.remove_response(inventory=inv, output="VEL", pre_filt=[0.005,0.01,45,50])
        sacfile = "%s/%s.%s.%s.%s.SAC" %(sac_path,net,sta,loc,cha)
        st.write(sacfile,format="SAC")
    except:
        print(net,sta)


def add_sac_header(sacfile,xmlfile_path,event_id,catalog,picks_file):
    os.putenv("SAC_DISPLAY_COPYRIGHT",'0')
    net = sacfile.split(".")[0]
    sta = sacfile.split(".")[1]
    loc = sacfile.split(".")[2]

    # stlo, stla
    xmlfile = xmlfile_path+"/%s.%s.xml" %(net,sta)
    inv = read_inventory(xmlfile).select(location=loc)
    stla = float(inv[0][0].latitude)
    stlo = float(inv[0][0].longitude)

    # evlo, evla, evdp
    command1 = "cat %s | awk '/%s/ {print}'" %(catalog, event_id)
    event_info = os.popen(command1).read().strip()
    evlo = float(event_info.split()[1])
    evla = float(event_info.split()[2])
    evdp = float(event_info.split()[3])

    # P(T0) and S(T1) wave arrival
    o = 0
    t0 = 0
    t1 = 0
    sta_match_string = "%s.%s" %(net,sta)
    command2 = "cat %s | awk '/%s/ {print}' | awk '/%s/ {print}'" \
               %(picks_file, event_id, sta_match_string)
    picks = os.popen(command2).read().strip().split("\n")
    if picks[0] != "":
        for pick in picks:
            pick = pick.strip()
            pick_sp = pick.split()
            phase = pick_sp[-1]
            origintime = pick_sp[1]
            o = datetime.datetime.strptime(origintime,"%Y-%m-%dT%H:%M:%S.%fZ")
            o_jday = o.strftime("%j")
            o_msec = int(o.microsecond/1000+0.5)
            if phase == "P":
                t0 = pick_sp[3]
                t0 = datetime.datetime.strptime(t0,"%Y-%m-%dT%H:%M:%S.%fZ")
                t0_jday = t0.strftime("%j")
                t0_msec = int(t0.microsecond/1000+0.5)
            elif phase == "S":
                t1 = pick_sp[3]
                t1 = datetime.datetime.strptime(t1,"%Y-%m-%dT%H:%M:%S.%fZ")
                t1_jday = t1.strftime("%j")
                t1_msec = int(t1.microsecond/1000+0.5)
    
    s = ""
    s += "readhdr %s.%s.*.SAC\n" %(net,sta)
    s += "ch stla %s\n" %(stla)
    s += "ch stlo %s\n" %(stlo)
    s += "ch evla %s\n" %(evla)
    s += "ch evlo %s\n" %(evlo)
    s += "ch evdp %s\n" %(evdp)
    s += "ch lcalda True\n" # calculate distance and azimuth
    if type(o) is datetime.datetime:
        s += "ch o gmt {} {} {} {} {} {}\n".format(o.year,o_jday,o.hour,
                                            o.minute,o.second,o_msec)
    if type(t0) is datetime.datetime:
        s += "ch t0 gmt {} {} {} {} {} {}\n".format(t0.year,t0_jday,t0.hour,
                                             t0.minute,t0.second,t0_msec)
    if type(t1) is datetime.datetime:
        s += "ch t1 gmt {} {} {} {} {} {}\n".format(t1.year,t1_jday,t1.hour,
                                             t1.minute,t1.second,t1_msec)
    s += "wh\n"
    s += "q\n"
    subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(s.encode())


main_path = os.getcwd()
xmlfile_path = main_path+"/stations"
miniseed_path = main_path+"/waveforms"
catalog = main_path+"/data/AACSE_catalog.dat"
picks_file = main_path+"/data/picks.dat"
dirs = os.listdir(miniseed_path)
print("===== converting miniseed to SAC =====")
for dirname in dirs:
    os.chdir("%s/%s" %(miniseed_path,dirname))
    print("processing on %s" %(dirname))
    sac_path = main_path+"/processedSeismograms/"+dirname
    mkdir(sac_path)
    miniseed_list = glob.glob("*.mseed")
    for miniseed in miniseed_list:
        mseedtosac(miniseed,xmlfile_path,sac_path)
print("===== Conversion completed =====")

print("===== Adding sac header =====")
seismograms = main_path+"/processedSeismograms"
dirs = os.listdir(seismograms)
for dirname in dirs:
    print("adding on %s" %(dirname))
    os.chdir("%s/%s" %(seismograms,dirname))
    sac_list = glob.glob("*.??Z.SAC")
    for sacfile in sac_list:
        add_sac_header(sacfile,xmlfile_path,dirname,catalog,picks_file)
print("===== completed =====")
    