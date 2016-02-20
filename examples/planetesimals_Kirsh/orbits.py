import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

diagnostics = 1

fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=',')
time = data[:,0]/6.28316

ms=3
if diagnostics == 1:
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,10))
    axes[0].plot(time,data[:,2], 'o', ms=ms, markeredgecolor='none')
    axes[0].set_ylabel('semi-major axis of Planet (AU, linear scale)')
    axes[1].plot(time,data[:,1], 'o', ms=ms, markeredgecolor='none')
    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    axes[1].set_xlim([1,max(data[:,0])])
    axes[1].set_ylabel('Energy (log scale)')
    axes[2].plot(data[:,0],data[:,3]-3, 'o', ms=ms, markeredgecolor='none')
    #axes[2].set_ylim([0,max(data[:,3]) + 10])
    axes[2].set_xlim([1,max(data[:,0])])
    axes[2].set_ylabel('N particles in global (log scale)')
    axes[2].set_xlabel('simulation time (yrs)')
    axes[2].set_xscale('log')
    axes[2].set_yscale('log')
    #axes[2].plot(data[:,0],data[:,5], 'o', ms=ms, markeredgecolor='none')
    #axes[2].set_xscale('log')
    #axes[2].set_yscale('log')
    #axes[2].set_ylabel('real life elapsed time (seconds)')
    #axes[2].set_xlabel('simulation time (yrs)')
    #axes[2].set_xlim([1,max(data[:,0])])
else:
    plt.plot(time,data[:,1], 'o', ms=ms, markeredgecolor='none')
    plt.plot(time,3e-10*data[:,0]**(0.5),color='black',label='t^1/2 growth')
    plt.ylabel('Energy')
    plt.xlabel('time (years)')
    plt.legend(loc='upper left',prop={'size':10})
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([0.5,data[-1,0]])

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
