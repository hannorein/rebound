#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random
#import numpy as np

#Specify what runs you want *****************************

#legend:[number of runs, N_planetesiamls, min arg1 val, max arg1 val, arg1 Type]
params = [50,50,0.1,16,'HSR']     #HSR
#params = [23,50,1e-5,0.5,'dt']    #dt

#*****************************

runs = []
log_output = (params[3]/params[2])**(1./params[0])
arg1 = params[2]
Np = str(params[1])
for i in xrange(0,params[0]):
    arg1 *= log_output
    seed = str(int(1000*random.random()))
    if params[4] == 'dt':
        arg1str = "{:.1e}".format(arg1)
    else:
        arg1str = "{:.2f}".format(arg1)
    name = 'output/'+params[4]+arg1str+'_Np'+Np+'_sd'+seed
    runs.append((arg1str,Np,seed,name))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+ ' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[runs[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
