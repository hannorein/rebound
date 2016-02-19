#This macro calculates the average energy for a set of PLANETESIMAL, SWIFTER and MERCURY runs and plots them all on one sexy figure.

import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import os

ms = 0.25
alpha = 0.5
choice = 1
names = ['time (years)','dE/E(0)','semi-major axis of planet(AU)']
outputn = ['time (years)','dE','a']

N_files = 0
dirP = str(sys.argv[1])
files = glob.glob(dirP+'*.txt')
data = []
n_it = 10e10

g_c = 0
b_c = 0
colors = []
colorsg = ['lime','chartreuse','forestgreen','springgreen','olivedrab','green']
colorsb = ['dodgerblue','blue','mediumblue','darkblue','royalblue','navy']

for f in files:
    try:
        split = f.split("_")
        if split[-2] == 'dt4.00':
            colors.append(colorsg[g_c])
            g_c += 1
        else:
            colors.append(colorsb[b_c])
            b_c += 1
        ff = open(f, 'r')
        lines = ff.readlines()
        length = len(lines)
        if length < n_it:   #need to find array with shortest length
            n_it = length
        data.append(lines)
        N_files += 1
    except:
        print 'couldnt read in data file '+f

E = np.zeros(shape=(N_files,n_it))
Eavg = np.zeros(n_it)
time = np.zeros(n_it)
vals_for_med = np.zeros(N_files)
for i in xrange(0,n_it):
    for j in range(0,N_files):
        split = data[j][i].split(",")
        vals_for_med[j] = float(split[choice])
        E[j][i] = vals_for_med[j]
    Eavg[i] = np.median(vals_for_med)
    time[i] = float(split[0])/6.283185

j=0
lc = len(colors)
for i in xrange(0,N_files):
    if colors[j] == colorsg[0]:
        plt.plot(time,E[i], ms=ms, color=colors[j], alpha=alpha, label = 'dt = 0.637')
    elif colors[j] == colorsb[0]:
        plt.plot(time,E[i], ms=ms, color=colors[j], alpha=alpha, label = 'dt = 2')
    else:
        plt.plot(time,E[i], ms=ms, color=colors[j], alpha=alpha)
    j += 1
    if j>=lc:
        j=0
plt.plot(time, Eavg, 'o', markeredgecolor='none', color='black', label='avg.')
if choice == 1:
    plt.plot(time,3e-10*time**(0.5),color='black')


##############################################
#Final plotting stuff
plt.legend(loc='lower right',prop={'size':10})
plt.ylabel(names[choice])
plt.xlabel('time (years)')
if choice == 1:
    plt.yscale('log')
    plt.xscale('log')
plt.xlim([0.5,time[-1]])
plt.title('Median Average '+names[choice]+' from '+str(N_files)+' files')
plt.savefig(dirP+'energy_avg_'+outputn[choice]+'.png')
plt.show()