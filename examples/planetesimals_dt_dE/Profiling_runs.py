#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
import time
#import numpy as np

#dt vs. dE for 1000 orbits - Np500
params=[(1e-5,500,12,"output/dt1e-5_Np500_sd12"),(2e-5,500,12,"output/dt2e-5_Np500_sd12"),(5e-5,500,12,"output/dt5e-5_Np500_sd12"),(1e-4,500,12,"output/dt1e-4_Np500_sd12"),(2e-4,500,12,"output/dt2e-4_Np500_sd12"),(5e-4,500,12,"output/dt5e-4_Np500_sd12"),(1e-3,500,12,"output/dt1e-3_Np500_sd12"),(2e-3,500,12,"output/dt2e-3_Np500_sd12"),(5e-3,500,12,"output/dt5e-3_Np500_sd12"),(1e-2,500,12,"output/dt1e-2_Np500_sd12"),(2e-3,500,12,"output/dt2e-3_Np500_sd12"),(5e-2,500,12,"output/dt5e-2_Np500_sd12"),(0.1,500,12,"output/dt0.1_Np500_sd12"),(0.2,500,12,"output/dt0.2_Np500_sd12"),(0.5,500,12,"output/dt0.5_Np500_sd12"),(1,500,12,"output/dt1_Np500_sd12"),(2,500,12,"output/dt2_Np500_sd12"),(5,500,12,"output/dt5_Np500_sd12"),(10,500,12,"output/dt10_Np500_sd12"),(20,500,12,"output/dt50_Np500_sd12"),(50,500,12,"output/dt50_Np500_sd12"),(100,500,12,"output/dt100_Np500_sd12"),(200,500,12,"output/dt200_Np500_sd12"),(500,500,12,"output/dt500_Np500_sd12")]

#dt vs. dE for 1000 orbits - Np50
params = [(1e-5,50,12,"output/dt1e-5_Np50_sd12"),(2e-5,50,12,"output/dt2e-5_Np50_sd12"),(5e-5,50,12,"output/dt5e-5_Np50_sd12"),(1e-4,50,12,"output/dt1e-4_Np50_sd12"),(2e-4,50,12,"output/dt2e-4_Np50_sd12"),(5e-4,50,12,"output/dt5e-4_Np50_sd12"),(1e-3,50,12,"output/dt1e-3_Np50_sd12"),(2e-3,50,12,"output/dt2e-3_Np50_sd12"),(5e-3,50,12,"output/dt5e-3_Np50_sd12"),(1e-2,50,12,"output/dt1e-2_Np50_sd12"),(2e-3,50,12,"output/dt2e-3_Np50_sd12"),(5e-2,50,12,"output/dt5e-2_Np50_sd12"),(0.1,50,12,"output/dt0.1_Np50_sd12"),(0.2,50,12,"output/dt0.2_Np50_sd12"),(0.5,50,12,"output/dt0.5_Np50_sd12"),(1,50,12,"output/dt1_Np50_sd12"),(2,50,12,"output/dt2_Np50_sd12"),(5,50,12,"output/dt5_Np50_sd12"),(10,50,12,"output/dt10_Np50_sd12"),(20,50,12,"output/dt50_Np50_sd12"),(50,50,12,"output/dt50_Np50_sd12"),(100,50,12,"output/dt100_Np50_sd12"),(200,50,12,"output/dt200_Np50_sd12"),(500,50,12,"output/dt500_Np50_sd12")]

length = len(params)

os.system('make')

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+ ' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[params[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
