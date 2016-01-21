#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
#import numpy as np

#params=[(10000000,50,13,"output/Np50_sd13.txt"),(10000000,50,26,"output/Np50_sd26.txt"),(10000000,50,34,"output/Np50_sd34.txt"),(10000000,50,51,"output/Np50_sd51.txt"),(10000000,50,63,"output/Np50_sd63.txt"),(10000000,50,65,"output/Np50_sd65.txt"),(10000000,50,70,"output/Np50_sd70.txt"),(10000000,50,10,"output/Np50_sd10.txt"),(5000000,100,16,"output/Np100_sd16.txt"),(5000000,100,18,"output/Np100_sd18.txt"),(5000000,100,28,"output/Np100_sd28.txt"),(5000000,100,37,"output/Np100_sd37.txt"),(5000000,100,52,"output/Np100_sd52.txt"),(5000000,100,62,"output/Np100_sd62.txt"),(5000000,100,67,"output/Np100_sd67.txt"),(5000000,100,72,"output/Np100_sd72.txt"),(1000000,500,11,"output/Np500_sd11.txt"),(1000000,500,12,"output/Np500_sd12.txt"),(1000000,500,22,"output/Np500_sd22.txt"),(1000000,500,33,"output/Np500_sd33.txt"),(1000000,500,49,"output/Np500_sd49.txt"),(1000000,500,62,"output/Np500_sd62.txt"),(1000000,500,65,"output/Np500_sd65.txt"),(1000000,500,71,"output/Np500_sd71.txt")]

params=[(5000000,100,16,"output/Np100_sd16.txt"),(5000000,100,18,"output/Np100_sd18.txt"),(5000000,100,28,"output/Np100_sd28.txt"),(5000000,100,37,"output/Np100_sd37.txt"),(5000000,100,52,"output/Np100_sd52.txt"),(5000000,100,62,"output/Np100_sd62.txt"),(5000000,100,67,"output/Np100_sd67.txt"),(5000000,100,72,"output/Np100_sd72.txt")]
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
