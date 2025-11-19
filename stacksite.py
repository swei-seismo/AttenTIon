import os
import glob
import math
import scipy.stats as stats
import matplotlib.pyplot as plt
from collections import defaultdict


def get_mean_ci(values):
    n = len(values)
    if n == 0:
        return None, None
    
    mean = sum(values)/n
    if n == 1:
        return mean, 0
    
    std_dev = math.sqrt(sum((x - mean)**2 for x in values)/(n-1))
    se = std_dev/math.sqrt(n)
    t_value = stats.t.ppf(0.975, df=n-1)
    ci = t_value*se
    
    return mean, ci


maindir = os.getcwd()
resultdir = maindir+"/workdir"
sitespecdir = maindir+"/data/sitespec"
plotsitespecdir = resultdir+"/plotsitespec"
for directory in [sitespecdir, plotsitespecdir]:
    os.makedirs(directory, exist_ok=True)

with open("stations.dat", "r") as f:
    for line in f:
        sta = line.strip().split()[0]
        file_pattern = "workdir/result/site/*_Presspec_%s.dat" %(sta)
        file_lst = glob.glob(file_pattern)
        if len(file_lst) <= 3:
            print(sta)
            continue
        else:
            data = defaultdict(lambda: {"residuals": [], "ratios": []})  
            for file in file_lst:
                with open(file, "r") as f1:
                    for line in f1:
                        parts = line.split()
                        freq, residual, ratio = float(parts[0]), float(parts[1]), float(parts[2])
                        data[freq]['residuals'].append(residual)
                        data[freq]['ratios'].append(ratio)

            results = []
            for freq in sorted(data.keys()):
                 avg_residual, ci_residual = get_mean_ci(data[freq]["residuals"])
                 avg_ratio, ci_ratio = get_mean_ci(data[freq]["ratios"])
                 results.append(f"{freq} {round(avg_residual, 4)} {round(ci_residual, 4)} {round(avg_ratio, 4)} {round(ci_ratio, 4)}")
            
            outputfile = sitespecdir+"/"+sta+"_Psitespec.dat"
            with open(outputfile, "w") as f2:
                f2.write("\n".join(results))
                   