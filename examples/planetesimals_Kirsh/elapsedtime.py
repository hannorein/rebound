#This macro calculates the average energy for a set of runs. Assumes that all the files are in energy_avg/.

import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import os
import re

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

twopi = 6.283185307

ms = 0.25
alpha = 0.5
colors = ['red','green','blue']

dirP = str(sys.argv[1])
files = glob.glob(dirP+'*_elapsedtime.txt')
files = sorted(files, key = natural_key)
time_array = ['4.00','8.00','12.76']
len_time = len(time_array)
HSR_array = ['0.25','0.50','1.00','2.00','4.00','8.00']
elapsedtime = np.zeros([len_time,len(HSR_array)])

for f in files:
    ff = open(f, 'r')
    lines = ff.readlines()
    split = lines[1]
    split = split.split()
    header = f.split("_")
    dt = header[-4]
    dt = re.sub('^dt', '', dt)
    HSR = header[-3]
    HSR = re.sub('^HSR', '', HSR)
    for i in xrange(0,len_time):
        for j in xrange(0,len(HSR_array)):
            if dt == time_array[i] and HSR == HSR_array[j]:
                elapsedtime[i][j] = float(split[-2])/3600

for i in xrange(0,len(HSR_array)):
    HSR_array[i] = float(HSR_array[i])

##############################################
#Final plotting stuff
for i in xrange(0,len_time):
    x = np.zeros(0)
    y = np.zeros(0)
    for j in xrange(0,len(HSR_array)):
        if elapsedtime[i][j] != 0:
            x = np.append(x,HSR_array[j])
            y = np.append(y,elapsedtime[i][j])
    plt.plot(x,y, ms = ms, color=colors[i], alpha=alpha, label = 'dt$_{sim}$='+str(round(100*float(time_array[i])/twopi)/100)+' years')
plt.legend(loc='upper left',prop={'size':10})
plt.ylabel('elapsed simulation time (hours)')
plt.xlabel('Hybrid Switch Ratio value')
plt.xlim([0,8.1])
plt.title('Elapsed time vs. HSR')
plt.savefig(dirP+'elapsed_time.png')
plt.show()