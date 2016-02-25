#This macro calculates the average energy for a set of PLANETESIMAL, SWIFTER and MERCURY runs and plots them all on one sexy figure.

import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import os

N_planetesimals = raw_input("Enter # planetesimals for runs you wish to average: ")
ms = 0.25
alpha = 0.5

planetesimal = 1
swifter = 1
mercury = 1

##############################################
#PLANETESIMAL
if planetesimal == 1:
    N_files = 0
    dirP = str(sys.argv[1])
    files = glob.glob(dirP+'*.txt')
    data = []
    n_it = 10e10
    for f in files:
        split = f.split("/")
        split2 = split[-1].split("_")
        if split2[0] == 'Np'+N_planetesimals:
            try:
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
            vals_for_med[j] = float(split[1])
            E[j][i] = vals_for_med[j]
        Eavg[i] = np.median(vals_for_med)
        time[i] = float(split[0])

    j=0
    colors = ['lime','chartreuse','forestgreen','springgreen','olivedrab','green']
    lc = len(colors)
    for i in xrange(0,N_files):
        plt.plot(time,E[i], ms=ms, color=colors[j], alpha=alpha)
        j += 1
        if j>=lc:
            j=0
    plt.plot(time, Eavg, 'o', markeredgecolor='none', color='Green', label='PLANETESIMAL avg.')
    plt.plot(time,3e-10*time**(0.5),color='black')

##############################################
#SWIFTER
if swifter == 1:
    print '...Finished Planetesimal, doing Swifter now...'
    dir = '../../../swifter/example/input_files/'
    files = [x[0] for x in os.walk(dir)]
    files = files[1:]
    N_files = 0
    n_it = 10e10
    print 'reading in data'
    data = []
    for f in files:
        split = f.split("/")
        split2 = split[-1].split("_")
        if split2[0] == 'Np'+N_planetesimals:
            try:
                ff = open(f+'/energyoutput.txt', 'r')
                lines = ff.readlines()
                length = len(lines)
                if length < n_it:   #need to find array with shortest length
                    n_it = length
                data.append(lines)
                N_files += 1
            except:
                print 'couldnt read in data file '+f+'/energyoutput.txt'

    E = np.zeros(shape=(N_files,n_it))
    Eavg = np.zeros(n_it)
    time = np.zeros(n_it)
    vals_for_med = np.zeros(N_files)
    print 'calculating avg energy'
    for i in xrange(1,n_it):
        for j in range(0,N_files):
            split = data[j][i].split()
            vals_for_med[j] = float(split[1])
            E[j][i] = vals_for_med[j]
        Eavg[i] = np.median(vals_for_med)
        time[i] = float(split[0])

    j=0
    colors = ['dodgerblue','blue','mediumblue','darkblue','royalblue','navy']
    lc = len(colors)
    for i in xrange(0,N_files):
        plt.plot(time,E[i], ms=ms, color=colors[j], alpha=alpha)
        j += 1
        if j>=lc:
            j=0
    plt.plot(time, Eavg, 'o', markeredgecolor='none', color='Blue', label='Swifter Avg.')

##############################################
#MERCURY
if mercury == 1:
    print '...Finished Swifter, doing Mercury now...'
    dir = '../../../mercury6/input_files/'
    files = [x[0] for x in os.walk(dir)]
    files = files[1:]
    N_files = 0
    n_it = 10e10

    print 'reading in data'
    data = []
    for f in files:
        split = f.split("/")
        split2 = split[-1].split("_")
        if split2[0] == 'Np'+N_planetesimals:
            try:
                ff = open(f+'/eo.txt', 'r')
                lines = ff.readlines()
                length = len(lines)
                if length < n_it:   #need to find array with shortest length
                    n_it = length
                data.append(lines)
                N_files += 1
            except:
                print 'couldnt read in data file '+f+'/eo.txt'

    E = np.zeros(shape=(N_files,n_it))
    Eavg = np.zeros(n_it)
    time = np.zeros(n_it)
    vals_for_med = np.zeros(N_files)

    print 'calculating avg energy'
    for i in xrange(0,n_it):
        split = data[0][i].split()
        vals_for_med[0] = float(split[1])
        E[0][i] = vals_for_med[0]
        t=float(split[0])
        for j in range(1,N_files):
            split = data[j][i].split()
            if float(split[0]) == t:
                vals_for_med[j] = float(split[1])
                E[j][i] = vals_for_med[j]
            else:
                print 'problem, times dont match'
        Eavg[i] = np.median(vals_for_med)
        time[i] = float(split[0])*0.01721420632

    j=0
    colors = ['lightcoral','firebrick','red','tomato','indianred','lightsalmon']
    lc = len(colors)
    for i in xrange(0,N_files):
        plt.plot(time,E[i], ms=ms, color=colors[j], alpha=alpha)
        j += 1
        if j>=lc:
            j=0
    plt.plot(time, Eavg, 'o', markeredgecolor='none', color='Red', label='MERCURY Avg.')

##############################################
#Final plotting stuff
plt.legend(loc='lower right',prop={'size':10})
plt.ylabel('dE/E(0)', fontsize=16)
plt.xlabel('time (years)', fontsize=16)
plt.yscale('log')
plt.xscale('log')
plt.xlim([0.5,time[-1]])
#plt.title('Median Average dE/E(0) from '+str(N_files)+' runs for Np'+N_planetesimals)
plt.savefig(dirP+'energy_avg_Np'+N_planetesimals+'_FULL.png')
plt.show()