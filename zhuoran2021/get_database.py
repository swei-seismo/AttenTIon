import os
import glob
import subprocess


def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def sactotxt(sacfile):
    os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
    net = sacfile.split(".")[0]
    sta = sacfile.split(".")[1]
    sac_list = glob.glob("%s.%s.*.SAC" %(net,sta))
    for sac in sac_list:
        txt_filename = sac.replace("SAC","txt")
        s = ""
        s += "r %s \n" %(sacfile)
        s += "rmean; rtr; taper\n"
        s += "ch allt (0 - &1,T0&) iztype IT0 \n" # set T0 as zero time
        s += "w over \n"
        s += "w alpha %s \n" %(txt_filename)
        s += "q \n"
        subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(s.encode())


main_path = os.getcwd()
eventname = main_path+"/data/Eventname"
dir_sacfiles = main_path+"/data/processedSeismograms"
dir_database = main_path+"/database"
mkdir(dir_database)

fevent = open(eventname,"r")
for event in fevent.readlines():
    event = event.strip()
    fdatabase = dir_database+"/database_%s" %(event)
    database = open(fdatabase,"w")
    database.write("net, sta, P_arrival, S_arrival, dist, baz, gcarc, delta, b \n")
    os.chdir(dir_sacfiles+"/"+event)
    sacfiles = glob.glob("*.??Z.SAC")
    for sacfile in sacfiles:
        command = "saclst knetwk kstnm T0 T1 dist baz gcarc \
                   delta b f %s" %(sacfile)
        sac_headers = os.popen(command).read().strip()
        net = sac_headers.split()[1]
        sta = sac_headers.split()[2]
        T0 = sac_headers.split()[3]
        T1 = sac_headers.split()[4]
        dist = sac_headers.split()[5]
        baz = sac_headers.split()[6]
        gcarc = sac_headers.split()[7]
        delta = sac_headers.split()[8]
        b = sac_headers.split()[9]
        database.write(net+" "+sta+" "+T0+" "+T1+" "+dist+" " \
                       +baz+" "+gcarc+" "+delta+" "+b+"\n")
        sactotxt(sacfile)
    database.close()
fevent.close()