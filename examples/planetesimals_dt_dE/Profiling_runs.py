#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
import time
#import numpy as np

#dt vs. dE for 1000 orbits - Np500
#params=[(1e-5,500,42,"output/dt1e-5_Np500_sd42"),(2e-5,500,42,"output/dt2e-5_Np500_sd42"),(5e-5,500,42,"output/dt5e-5_Np500_sd42"),(1e-4,500,42,"output/dt1e-4_Np500_sd42"),(2e-4,500,42,"output/dt2e-4_Np500_sd42"),(5e-4,500,42,"output/dt5e-4_Np500_sd42"),(1e-3,500,42,"output/dt1e-3_Np500_sd42"),(2e-3,500,42,"output/dt2e-3_Np500_sd42"),(5e-3,500,42,"output/dt5e-3_Np500_sd42"),(1e-2,500,42,"output/dt1e-2_Np500_sd42"),(2e-2,500,42,"output/dt2e-2_Np500_sd42"),(5e-2,500,42,"output/dt5e-2_Np500_sd42"),(0.1,500,42,"output/dt0.1_Np500_sd42"),(0.2,500,42,"output/dt0.2_Np500_sd42"),(0.5,500,42,"output/dt0.5_Np500_sd42")]

#dt vs. dE for 1000 orbits - Np50
#params = [(1e-5,50,42,"output/dt1e-5_Np50_sd42"),(2e-5,50,42,"output/dt2e-5_Np50_sd42"),(5e-5,50,42,"output/dt5e-5_Np50_sd42"),(1e-4,50,42,"output/dt1e-4_Np50_sd42"),(2e-4,50,42,"output/dt2e-4_Np50_sd42"),(5e-4,50,42,"output/dt5e-4_Np50_sd42"),(1e-3,50,42,"output/dt1e-3_Np50_sd42"),(2e-3,50,42,"output/dt2e-3_Np50_sd42"),(5e-3,50,42,"output/dt5e-3_Np50_sd42"),(1e-2,50,42,"output/dt1e-2_Np50_sd42"),(2e-2,50,42,"output/dt2e-2_Np50_sd42"),(5e-2,50,42,"output/dt5e-2_Np50_sd42"),(0.1,50,42,"output/dt0.1_Np50_sd42"),(0.2,50,42,"output/dt0.2_Np50_sd42"),(0.5,50,42,"output/dt0.5_Np50_sd42"),(1,50,42,"output/dt1_Np50_sd42")]

#HSR vs. dE for 1Myr orbits - Np50 - HSR = argument 1!
params = [(0.1,50,22,"output/HSR0.1_Np50_sd22"),(0.25,50,22,"output/HSR0.25_Np50_sd22"),(0.5,50,22,"output/HSR0.5_Np50_sd22"),(0.75,50,22,"output/HSR0.75_Np50_sd22"),(1,50,22,"output/HSR1_Np50_sd22"),(1.5,50,22,"output/HSR1.5_Np50_sd22"),(2,50,22,"output/HSR2_Np50_sd22"),(2.5,50,22,"output/HSR2.5_Np50_sd22"),(3,50,22,"output/HSR3_Np50_sd22"),(4,50,22,"output/HSR4_Np50_sd22"),(6,50,22,"output/HSR6_Np50_sd22"),(8,50,22,"output/HSR8_Np50_sd22"),(12,50,22,"output/HSR12_Np50_sd22"),(16,50,22,"output/HSR16_Np50_sd22")]

#HSR vs. dE for 1Myr orbits - Np500 - HSR = argument 1!
#params = [(0.1,500,22,"output/HSR0.1_Np500_sd22"),(0.25,500,22,"output/HSR0.25_Np500_sd22"),(0.5,500,22,"output/HSR0.5_Np500_sd22"),(0.75,500,22,"output/HSR0.75_Np500_sd22"),(1,500,22,"output/HSR1_Np500_sd22"),(1.5,500,22,"output/HSR1.5_Np500_sd22"),(2,500,22,"output/HSR2_Np500_sd22"),(2.5,500,22,"output/HSR2.5_Np500_sd22"),(3,500,22,"output/HSR3_Np500_sd22"),(4,500,22,"output/HSR4_Np500_sd22"),(6,500,22,"output/HSR6_Np500_sd22"),(8,500,22,"output/HSR8_Np500_sd22"),(12,500,22,"output/HSR12_Np500_sd22"),(16,500,22,"output/HSR16_Np500_sd22")]

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
