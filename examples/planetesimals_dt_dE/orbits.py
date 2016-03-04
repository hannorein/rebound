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

ms=3
if diagnostics == 1:
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(10,10))
    axes[0].plot(data[:,0],data[:,1], 'o', ms=ms, markeredgecolor='none')
    axes[0].plot(data[:,0],3e-11*data[:,0]**(0.5),color='black',label='t^1/2 growth')
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    #axes[0].set_xlim([20,30])
    axes[0].set_xlim([1,max(data[:,0])])
    axes[0].set_ylabel('Energy')
    axes[1].plot(data[:,0],data[:,2], 'o', ms=ms, markeredgecolor='none', label='global: Np')
    axes[1].plot(data[:,0],data[:,3], 'o', ms=ms, markeredgecolor='none', label='mini: Np')
    axes[1].legend(loc='upper left',prop={'size':10})
    axes[1].set_yscale('log')
    #axes[1].set_ylim([0,max(data[:,2])])
    #axes[1].set_ylim([0,4])
    axes[2].plot(data[:,0],data[:,4], 'o', ms=ms, markeredgecolor='none')
    axes[2].set_yscale('log')
    axes[2].set_ylabel('real life elapsed time (seconds)')
    axes[2].set_xlabel('simulation time (yrs)')
else:
    plt.plot(data[:,0],data[:,1], 'o', ms=ms, markeredgecolor='none')
    plt.plot(data[:,0],3e-10*data[:,0]**(0.5),color='black',label='t^1/2 growth')
    plt.ylabel('Energy')
    plt.xlabel('time (years)')
    plt.legend(loc='upper left',prop={'size':10})
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([0.5,data[-1,0]])

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
