#!/usr/local/opt/software/python/python-3.5.1-gnu/bin/python3

# Path to python3 installation on the USQHPC

# to find this path for a HPC installation:
#	module avail pyth
#	module load python/3.5.1-gnu		# or whatever version needed
#	which python 						# gives you the path
#------------------------------------------------------------------------------
# Change internal python path to local user directory for python modules. 
# Should have the same version as the python installation in the shebang.

import sys
sys.path.insert(0,'/home/u1092862/.local/lib/python3.5/site-packages')

# install new packages with comand for local user (example for rebound):
# pip install --proxy proxy.usq.edu.au:8080 --upgrade --user rebound
