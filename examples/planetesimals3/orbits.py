import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=',')
plt.plot(data[:,0],data[:,1], 'o', ms=msval, markeredgecolor='none')
plt.plot(data[:,0],3e-10*data[:,0]**(0.5),color='black',label='t^1/2 growth')

plt.legend(loc='upper left',prop={'size':10})
plt.xscale('log')

plt.ylabel(names[arg1])
plt.xlabel('time (years)')
plt.yscale('log')
plt.xlim([0.5,data[-1,0]])
plt.show()
