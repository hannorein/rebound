#This macro makes a plot of dt vs. energy (and also dt vs. time)

import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import os
import re

def sort(x, y):
    length = len(x)
    for i in xrange(0,length):
        for j in xrange(i,length):
            if x[j] < x[i]:
                tempx = x[j]
                tempy0 = y[0][j]
                tempy1 = y[1][j]
                x[j] = x[i]
                x[i] = tempx
                y[0][j] = y[0][i]
                y[0][i] = tempy0
                y[1][j] = y[1][i]
                y[1][i] = tempy1
    return x, y

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

#plot dt vs.
y_choice = int(sys.argv[1])      #0 = plot elapsed time, 1 = energy
x_choice = 0                     #0 = dt, 1 = HSR

dirP = str('dtvdE_files/')
Navg = 10                       #number of points at end of .txt file to average energy over

ext = 'dt'
xname = 'timestep (years)'
xvals = [1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3,2e-3,5e-3,0.01,0.05,0.1,0.2,0.5]
yvals = ['Np50_sd12','Np50_sd11','Np50_sd42','Np500_sd12','Np500_sd11','Np500_sd42']
marker = ['^--', '^--','^--','o-', 'o-','o-']
label = ['Np50','','','Np500','','']
if x_choice == 1:
    ext = 'HSR'
    xname = 'HSR (hill radii)'
    xvals = [0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,4,6,8,12,16]

lenx = len(xvals)
leny = len(yvals)
for i in xrange(0,leny):
    files = glob.glob(dirP+'*'+yvals[i]+'_elapsedtime.txt')
    files = sorted(files, key = natural_key)
    ET = np.zeros(lenx)
    dE = np.zeros(lenx)
    for f in files:
        header = f.split("_")
        xx = header[-4]
        x = float(re.sub('^files/'+ext, '', xx))
        for j in xrange(0,lenx):
            if x == xvals[j]:
                index = j
        ff = open(f, 'r')
        lines = ff.readlines()
        elapsed = lines[1]
        elapsed = elapsed.split()
        ET[index] = float(elapsed[-2])/3600
        txtfile = f.split("_elapsedtime.txt")
        fff = open(txtfile[0]+".txt","r")
        lines = fff.readlines()
        length = len(lines)
        Emed = np.zeros(Navg)
        for k in xrange(0,Navg):
            split = lines[-1-k]
            split = split.split(",")
            Emed[k] = float(split[1])
        med = np.median(Emed)
        if med != med:
            med = 1
        dE[index] = med
    if y_choice == 0:
        y = ET
    else:
        y = dE
    plt.plot(xvals, y, marker[i],label=label[i])

if y_choice == 0:
    name = 'elapsed time (seconds)'
    oname = 'ET'
else:
    name = 'dE/E(0)'
    oname = 'dE'
#x, y = sort(x, y)
plt.yscale('log')
if x_choice == 0:
    plt.xscale('log')
plt.ylabel(name)
plt.xlabel(xname)
plt.title('Integrating 2p, 2pl system for 1000 orbits')
#plt.ylim([np.min(y),2])
plt.legend(loc='lower right',prop={'size':10})
plt.savefig('dtvdE_files/'+ext+'_v_'+oname+'.png')
plt.show()

