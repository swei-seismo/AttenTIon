import os

def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


main_path = os.getcwd()
seismograms_path = main_path+"/data/processedSeismograms"
mkdir(seismograms_path)

backup_seis_path = main_path+"/processedSeismograms"
dirs = os.listdir(backup_seis_path)
for dirname in dirs:
    os.chdir("%s/%s" %(backup_seis_path,dirname))
    mkdir("%s/%s" %(seismograms_path,dirname))
    command = "saclst t0 t1 f *.SAC | awk -F ' ' ' \
               {if($2!=-12345) print $1}'"
    sacfiles = os.popen(command).read().strip().split("\n")
    for sacfile in sacfiles:
        os.system("cp ./%s %s/%s" %(sacfile,seismograms_path,dirname))