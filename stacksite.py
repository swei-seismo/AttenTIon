#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Calculate the averaged site effects for each station
And plot the individual site effect with the averaged site effect for each station
Need to input 
    1) df and endf values
    2) resultdir
    3) stafl: station name at the first column

...note:
    curretnly only P site spectrum is considered and calculated 
'''
import os, math
from matplotlib import pyplot as plt
import numpy as np
import glob
df = 0.0391
endf = 20

maindir = os.getcwd()
resultdir = '/mnt/scratch/yurong/tstar/workdir_nosite/result/'
sitespecdir = maindir + '/data/sitespec/'
plotsitespecdir = sitespecdir + '/plotsitespec/'
if not os.path.isdir(plotsitespecdir):
    os.makedirs(plotsitespecdir)

stalst = []


stafl = "./data/stations.lst"
for line in open(stafl).readlines():
    stalst.append(line.split()[0])

fre = np.linspace(df,endf,math.ceil(endf/df))
for i in range(len(fre)):
    # truncate the frequency numbers to 4 decimal places without rounding
    fre[i]=int(fre[i]*10000)/10000
    #rounding:
    #fre[i]=format(fre[i],'.4f')
    
for sta in stalst:
    allspec = [0.00] * math.ceil(endf/df)
    num = [0] * math.ceil(endf/df)
    # x, y save frequencies and spectra for each record
    x = []
    y = []
    print(sta)
    resfl = resultdir+'*%s.dat'%sta
    reslst = glob.glob(resfl)
    if len(reslst) == 0:
        print("no site spectra for %s"%sta)
        continue
    if len(reslst) <= 8:
        print("less then 8 site spectra for %s"%sta)
        continue
    print("%d residual spectra for %s"%(len(reslst),sta))
    for i,sitefl in enumerate(reslst):
        ifre = []
        ispec = []
        for line in open(sitefl).readlines():
            f = float(line.split()[0])
            allspec[(round((f/df)-1))] += float(line.split()[1])
            num[(round((f/df)-1))] += 1
            ifre.append(float(line.split()[0]))
            ispec.append(float(line.split()[1]))
        x.append(ifre)
        y.append(ispec)

    avespec = list(np.array(allspec)/np.array(num))

    outfl = open(sitespecdir + 'Psite_%s.spec'%sta, 'w')
    for i in range(len(fre)):
        outfl.write("%f %f \n"%(fre[i], avespec[i]))
    outfl.close()

    plt.figure()
    plt.clf()

    for i in range(len(x)):
        plt.plot(x[i],y[i],'lightgray')
        
    plt.plot(fre,avespec,c='red', linewidth=2) 
    plt.title('P wave Site Effect: %s %d traces'%(sta,len(x)),fontsize=14)
    plt.xlabel('Frequency(Hz)',fontsize=11)
    plt.ylabel('Residual Spectrum(nm/s)',fontsize=11)
    plt.savefig(plotsitespecdir + sta+'_site.pdf')

    


    
        
