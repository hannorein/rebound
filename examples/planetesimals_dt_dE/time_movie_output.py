import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import math
pi = math.pi
from mpl_toolkits.mplot3d import Axes3D
import re

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

dir = 'xyz_outputs/'
files = glob.glob(dir+'*.txt')
files = sorted(files, key=natural_key)
N_files = len(files)  #number of files we're dealing with

#colors - need to improve this later for any number of active/passive bodies
N_planets = 500
colors =  ["black" for x in range(N_planets)]
colors[0] = 'yellow'
colors[1] = 'orange'
colors[2] = 'red'

for i in xrange(0,N_files):
    fos = open(files[i], 'r')
    t,id,x,y,z,dE = np.loadtxt(fos, delimiter=',', unpack='True')
    dx = x[1] - x[46]
    dy = y[1] - y[46]
    dz = z[1] - z[46]
    r12 = (dx*dx + dy*dy + dz*dz)**0.5
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z,c=colors[0:len(x)], lw=0)
    ax.scatter(x[0:3],y[0:3],z[0:3],c=colors[0:3],lw=0, s=30)
    ax.set_title('t='+str(t[0])+', r_0,46='+str(r12)+', dE='+str(dE[0]))
    ax.view_init(elev = 90, azim=100)    #both are in degrees. elev = 0 or 90 is what you want
    output = files[i]
    output = re.sub('\.txt$', '', output)
    plt.savefig(output+'.png')
    print 'completed iteration',i