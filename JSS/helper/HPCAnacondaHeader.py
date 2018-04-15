#!/home/u1092862/.local/anaconda3/envs/rebound/bin/python

# Path to local user anaconda python3 on USQHPC

# to find this path 
#	conda env list						# active environment is marked with *
#	source activate rebound				# activate needed environment
#	which python 						# gives you the path
#------------------------------------------------------------------------------
# Change internal python path to local user directory for python modules. 
# Should have the same version as the python installation in the shebang.

import sys
sys.path.insert(0,'/home/u1092862/.local/anaconda3/envs/rebound/lib/python3.6/site-packages')

# install new packages with comand for local anaconda (example for rebound):
# 	pip install --proxy proxy.usq.edu.au:8080 --upgrade rebound
