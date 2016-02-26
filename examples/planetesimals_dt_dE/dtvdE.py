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
plot_choice = int(sys.argv[1])      #0 = plot elapsed time, 1 = energy

N_files = 0
dirP = str('dtvdE_files/')
files = glob.glob(dirP+'*Np500_*_elapsedtime.txt')
files = sorted(files, key = natural_key)
size = (2,len(files))
ET = np.zeros(size)
dE = np.zeros(size)
dt = np.zeros(len(files))
Navg = 10

i=0
for f in files:
    header = f.split("_")
    dtt = header[-4]
    dt[i] = float(re.sub('^files/dt', '', dtt))
    temp = header[-3]
    filenames = [f, re.sub('Np500','Np50',f)]
    for index in xrange(0,2):
        try:
            ff = open(filenames[index], 'r')
            lines = ff.readlines()
            elapsed = lines[1]
            elapsed = elapsed.split()
            ET[index][i] = float(elapsed[-2])/3600
            txtfile = filenames[index].split("_elapsedtime.txt")
            fff = open(txtfile[0]+".txt","r")
            lines = fff.readlines()
            length = len(lines)
            Emed = np.zeros(Navg)
            for j in xrange(0,Navg):
                split = lines[-1-j]
                split = split.split(",")
                Emed[j] = float(split[1])
            med = np.median(Emed)
            if med != med:
                med = 1
            dE[index][i] = med
        except:
            print 'couldnt find file: '+filenames[index]
    i+=1

if plot_choice == 0:
    y = ET
    name = 'elapsed time (seconds)'
    oname = 'ET'
else:
    y = dE
    name = 'dE/E(0)'
    oname = 'dE'
dt, y = sort(dt, y)
plt.plot(dt, y[0], 'o-',label='Np = 500')
plt.plot(dt, y[1], 'o-',label='Np = 50')
plt.yscale('log')
plt.xscale('log')
plt.ylabel(name)
plt.xlabel('timestep (years)')
plt.title('Integrating 2p, 2pl system for 1000 orbits')
#plt.ylim([np.min(y),2])
plt.legend(loc='lower right',prop={'size':10})
plt.savefig('dtvdE_files/dt_v_'+oname+'.png')
plt.show()

