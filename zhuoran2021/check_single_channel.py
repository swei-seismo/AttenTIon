import os
import glob

main_path = os.getcwd()
backup_seis_path = main_path+"/processedSeismograms"
dirs = os.listdir(backup_seis_path)
for dirname in dirs:
    os.chdir("%s/%s" %(backup_seis_path,dirname))
    sacfiles = glob.glob("*.??Z.SAC")
    for sacfile in sacfiles:
        net = sacfile.split(".")[0]
        sta = sacfile.split(".")[1]
        lst = glob.glob("%s.%s.*.SAC" %(net,sta))
        if len(lst) == 3:
            continue
        else:
            loc = sacfile.split(".")[2]
            chan = sacfile.split(".")[3]
            horizontal_one = chan.replace("Z","E")
            horizontal_two = chan.replace("Z","N")
            horizontal_sacfile_one = "%s.%s.%s.%s.SAC" %(net,sta,loc,horizontal_one)
            horizontal_sacfile_two = "%s.%s.%s.%s.SAC" %(net,sta,loc,horizontal_two)
            os.system("cp %s %s" %(sacfile,horizontal_sacfile_one))
            os.system("cp %s %s" %(sacfile,horizontal_sacfile_two))
