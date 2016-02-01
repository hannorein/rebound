#This macro calculates the difference in orbital growth between integrators.

import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import re

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
        RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

dirhybarid = 'xyz_outputs/'
dirswifter = '../../../swifter/example/xyz_outputs/'
filesh = glob.glob(dirhybarid+'*.txt')
filesh = sorted(filesh, key = natural_key)
filess = glob.glob(dirswifter+'*.txt')
filess = sorted(filess, key = natural_key)

fh = open(filesh[0],'r')
t,id,xh,yh,zh = np.loadtxt(fh, delimiter=',', unpack='True')
Nbods = len(xh)
N_files = len(filess)
if len(filesh) < N_files:
    N_files = len(filesh)

dr = np.zeros((N_files,Nbods + 1))
time = np.zeros(N_files)

for i in xrange(0,N_files):
    fh = open(filesh[i],'r')
    t,id,xh,yh,zh = np.loadtxt(fh, delimiter=',', unpack='True')
    fs = open(filess[i],'r')
    t,id,xs,ys,zs = np.loadtxt(fs, unpack='True')
    Nbods = len(xh)
    for j in xrange(0,Nbods):
        dx = xh[j] - xs[j]
        dy = yh[j] - ys[j]
        dz = zh[j] - zs[j]
        dr[i,j+1] = (dx*dx + dy*dy + dz*dz)**0.5
        dr[i,0] += dr[i,j+1]  #total error growth
    time[i] = t[0]

cmap = get_cmap(Nbods)
for i in xrange(0,Nbods):
    plt.plot(time,dr[:,i+1], ms=0.5, color=cmap(i), alpha=0.75, label='body '+str(i))
plt.plot(time, dr[:,0], 'o', markeredgecolor='none', color='black', label='Total Error')
plt.legend(loc='upper left',prop={'size':10})
plt.ylabel('dr deviation of swifter from hybarid')
plt.xlabel('time (years)')
plt.yscale('log')
plt.xscale('log')
#plt.xlim([0.5,time[-1]])
plt.title('Integrator deviations between Hybarid/Swifter for '+str(Nbods)+' bodies')
plt.savefig(dirhybarid+'integrator_drcomp_N'+str(Nbods)+'.png')
plt.show()