# -*- coding: utf-8 -*-

"""
This file contains routines to interface with MERCURY and SWIFTER.
It is not recommended to use these functions for anything other than testing.

Hanno Rein, Daniel Tamayo
2015

"""
import os
import tempfile
import shutil
import time
import math
import rebound

integrator_package = "REBOUND"
integrator_fullname = ""    

tmpdir = None
def reset_debug():
    global tmpdir
    if tmpdir:
        shutil.rmtree(tmpdir)
        tmpdir = None

def integrate_other_package(tmax,exactFinishTime=1):
    global tmpdir
    if integrator_package == "SWIFTER":
        _particles = rebound.module.particles
        oldwd = os.getcwd()
        paramin = """!
! Parameter file for the CHO run of the 4 giant planets and Pluto.
!
!NPLMAX         -1                 ! not used
!NTPMAX         -1                 ! not used
T0             %.16e
TSTOP          %.16e               ! length of simulation 
DT             %.16e               ! stepsize
PL_IN          pl.in
TP_IN          tp.in
IN_TYPE        ASCII
ISTEP_OUT      %d                   ! output every 10K years
BIN_OUT        bin.dat
OUT_TYPE       REAL8                ! 4-byte XDR formatted output
OUT_FORM       XV                  ! cartesian output 
OUT_STAT       NEW
ISTEP_DUMP     1000000000        ! system dump also every 10K years
J2             0.0E0               ! no J2 term
J4             0.0E0               ! no J4 term
CHK_CLOSE      no                  ! don't check for planetary close encounters
CHK_RMIN       -1.0                ! don't check for close solar encounters
CHK_RMAX       1000.0              ! discard outside of 1000 AU
CHK_EJECT      -1.0                ! ignore this check
CHK_QMIN       -1.0                ! ignore this check
!CHK_QMIN_COORD HELIO               ! commented out here
!CHK_QMIN_RANGE 1.0 1000.0          ! commented out here
ENC_OUT        enc.dat
EXTRA_FORCE    no                  ! no extra user-defined forces
BIG_DISCARD    yes                 ! output all planets if anything discarded
RHILL_PRESENT  no                  ! no Hill's sphere radii in input file
""" % ( rebound.module.t, tmax, rebound.module.dt,max(1,int((tmax-rebound.module.t)/rebound.module.dt)))

        if not tmpdir:
        # first call
            tmpdir = tempfile.mkdtemp()
            for f in ["swifter_whm", "swifter_tu4","swifter_symba","swifter_helio"]:
                os.symlink(oldwd+"/../../others/swifter/bin/"+f,tmpdir+"/"+f)
            os.symlink(oldwd+"/../../others/swifter/bin/tool_follow",tmpdir+"/tool_follow")
        os.chdir(tmpdir)
        smallin = """0\n"""
        with open("tp.in", "w") as f:
            f.write(smallin)
        bigin = """ %d
""" %(rebound.module.N)
        _G = rebound.module.G 
        for i in range(0,rebound.module.N):
            bigin += """%d    %.18e 
%.18e %.18e %.18e
%.18e %.18e %.18e
""" %(i+1, _G*_particles[i].m,_particles[i].x, _particles[i].y, _particles[i].z, _particles[i].vx, _particles[i].vy, _particles[i].vz)
        with open("pl.in", "w") as f:
            f.write(bigin)
        with open("param.in", "w") as f:
            f.write(paramin)
        #with open("param.dmp", "w") as f:
        #    f.write(paramin)
        starttime = time.time()    
        if integrator_fullname.lower() == "swifter-whm":
            os.system("echo param.in | ./swifter_whm > /dev/null")
        elif integrator_fullname.lower() == "swifter-symba":
            os.system("echo param.in | ./swifter_symba > /dev/null")
        elif integrator_fullname.lower() == "swifter-helio":
            os.system("echo param.in | ./swifter_helio > /dev/null")
        elif integrator_fullname.lower() == "swifter-tu4":
            os.system("echo param.in | ./swifter_tu4 > /dev/null")
        else:
            print("Integrator not found! %s"%integrator_fullname)
        endtime = time.time()    
        c_double.in_dll(rebound.module.clibrebound,"timing").value = endtime-starttime
    
        for i in range(1,rebound.module.N):
            with open("outputparams.txt", "w") as f:
                f.write("dump_param1.dat\n")
                f.write("{0}\n".format(i+1))
                f.write("1\n")
            try:
                os.system("./tool_follow < outputparams.txt > /dev/null")
                with open("follow.out", "r") as f:
                    lines = f.readlines()
                    line = lines[-1].split()
                    t = float(line[0].strip())
                    rebound.module.t = t
                    _particles[i].x = float(line[2].strip())
                    _particles[i].y = float(line[3].strip())
                    _particles[i].z = float(line[4].strip())
                    _particles[i].vx = float(line[5].strip())
                    _particles[i].vy = float(line[6].strip())
                    _particles[i].vz = float(line[7].strip())
            except:
                print("Something went wrong. Ignoring it for now. (%s)"%integrator_fullname)
                pass
        os.system("rm bin.dat")
        os.chdir(oldwd)
    if integrator_package == "MERCURY":
        k = 0.01720209895    
        facTime = math.sqrt(rebound.module.G)/k
        _particles = rebound.module.particles
        oldwd = os.getcwd()
        paramin = """)O+_06 Integration parameters  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------------------------
) Important integration parameters:
)---------------------------------------------------------------------
 algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = mvs
 start time (days)= 0.
 stop time (days) = %.16e
 output interval (days) = 1000000.1d0
 timestep (days) = %.16e
 accuracy parameter=1.d-12
)---------------------------------------------------------------------
) Integration options:
)---------------------------------------------------------------------
 stop integration after a close encounter = no
 allow collisions to occur = no
 include collisional fragmentation = no
 express time in days or years = days
 express time relative to integration start time = yes
 output precision = high
 < not used at present >
 include relativity in integration= no
 include user-defined force = no
)---------------------------------------------------------------------
) These parameters do not need to be adjusted often:
)---------------------------------------------------------------------
 ejection distance (AU)= 100
 radius of central body (AU) = 0.005
 central mass (solar) = %.16e
 central J2 = 0
 central J4 = 0
 central J6 = 0
 < not used at present >
 < not used at present >
 Hybrid integrator changeover (Hill radii) = 3.
 number of timesteps between data dumps = 50000000
 number of timesteps between periodic effects = 100000000
""" % ( tmax*facTime, rebound.module.dt*facTime,_particles[0].m)
        if not tmpdir:
            # first call
            tmpdir = tempfile.mkdtemp()
            for f in ["mercury", "message.in","files.in"]:
                os.symlink(oldwd+"/../../others/mercury6/"+f,tmpdir+"/"+f)
            os.chdir(tmpdir)
            smallin = """)O+_06 Small-body initial data  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------------------------
 style (Cartesian, Asteroidal, Cometary) = Ast
)---------------------------------------------------------------------
"""
            with open("small.in", "w") as f:
                f.write(smallin)
            bigin = """)O+_06 Big-body initial data  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------------------------
 style (Cartesian, Asteroidal, Cometary) = Cartesian
 epoch (in days) = %.16e
)---------------------------------------------------------------------
""" %(rebound.module.t*facTime)
            for i in range(1,rebound.module.N):
                bigin += """PART%03d    m=%.18e r=20.D0 d=5.43
 %.18e %.18e %.18e
 %.18e %.18e %.18e
  0. 0. 0.
""" %(i, _particles[i].m, _particles[i].x, _particles[i].y, _particles[i].z, _particles[i].vx/facTime, _particles[i].vy/facTime, _particles[i].vz/facTime)
            with open("big.in", "w") as f:
                f.write(bigin)
            with open("param.in", "w") as f:
                f.write(paramin)
        else:
            # Not first call
            os.chdir(tmpdir)
            with open("param.dmp", "w") as f:
                f.write(paramin)
        starttime = time.time()    
        #os.system("./mercury ")
        os.system("./mercury >/dev/null")
        endtime = time.time()    
        c_double.in_dll(rebound.module.clibrebound,"timing").value = endtime-starttime
        with open("big.dmp", "r") as f:
            lines = f.readlines()
            t= float(lines[4].split("=")[1].strip())
            rebound.module.t = t/facTime
            j = 1
            for i in range(6,len(lines),4):
                pos = lines[i+1].split()
                _particles[j].x = float(pos[0])
                _particles[j].y = float(pos[1])
                _particles[j].z = float(pos[2])
                vel = lines[i+2].split()
                _particles[j].vx = float(vel[0])*facTime
                _particles[j].vy = float(vel[1])*facTime
                _particles[j].vz = float(vel[2])*facTime
                j += 1

        os.chdir(oldwd)


# Import at the end to avoid circular dependence
from .simulation import Simulation
