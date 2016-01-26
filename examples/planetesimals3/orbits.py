import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=',')
plt.plot(data[:,0],data[:,1], 'o', ms=2, markeredgecolor='none')
plt.plot(data[:,0],3e-10*data[:,0]**(0.5),color='black',label='t^1/2 growth')

plt.legend(loc='upper left',prop={'size':10})
plt.xscale('log')

plt.ylabel('Energy')
plt.xlabel('time (years)')
plt.yscale('log')
plt.xlim([0.5,data[-1,0]])
file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
