#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
import time
#import numpy as np

#params=[(4,10,0.25),(4,10,0.5),(4,10,1),(4,10,2),(4,10,4),(4,10,8),(4,10,16),(4,20,2),(4,30,2),(4,40,2),(4,50,2),(4,60,2),(12.76,10,0.25),(12.76,10,0.5),(12.76,10,1),(12.76,10,2),(12.76,10,4),(12.76,10,8),(12.76,10,16),(12.76,20,2),(12.76,30,2),(12.76,40,2),(12.76,50,2),(12.76,60,2)]
params=[(4,10,0.25),(4,10,0.5),(4,10,1),(4,10,2),(4,10,4),(4,10,8),(8,11,0.25),(8,11,0.5),(8,11,1),(8,11,2),(8,11,4),(8,11,8),(12.76,12,0.25),(12.76,12,0.5),(12.76,12,1),(12.76,12,2),(12.76,12,4),(12.76,12,8)]

length = len(params)

os.system('make')

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[params[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
