import rebound as rb
import addons
from argparse import ArgumentParser as ap
import configparser
from astropy.time import Time
import datetime
import sys
from pathlib import Path

# parse commandline arguments
parser = ap(description='Reads command line arguments and builds configuration '\
			'file. With no other arguments provided, the configuration '	\
			'file for <project_name> is used. If there is no configuration '	\
			'file in the data directory, the program will create a new config '	\
			'file using the defaults as well as arguments if provided. -n will '\
		    'overwrite an existing file without warning otherwise an existing '	\
		    'configuration file is updated with the arguments provided' )
parser.add_argument('project_name', 
					type		= str,
					help		= 'Id of the project. (example: JSS-LR-12S for '\
								  '"Jupiter-Satellite-System-Long-Run-12-Satell'\
								  'ites). <project_name> is used in all project'\
								  ' related files.  Use whatever rocks your '	\
								  'boat! ')
parser.add_argument('-n', '--new', 											
					action		= 'store_const',
					const 		= True,
					help 		= 'writes a new ini file with the default '		\
								  'settings')
parser.add_argument('-r', '--rebuild', 											
					action		= 'store_const',
					const 		= True,
					help 		= 'force rebuild binary file. Use this one '	\
								  'after you have changed the config file by'	\
								  'editing outside this program')
parser.add_argument('-l', '--list_of_bodies', 
	                action		= 'store', 
	                type 		= str, 
	                nargs 		= '*', 
					default 	= ['Jupiter', 'Metis', 'Adrastea', 'Amalthea',
					              'Thebe', 'Io', 'Europa', 'Ganymede', 'Callisto',
					              'Sun', 'Saturn barycenter', 'Uranus barycenter', 
					              'Neptune barycenter'],
					help		= 'List of bodies for the model, separated by '	\
								  'blanks.')
parser.add_argument('-a', '--add_to_sun', 
	                action		= 'store', 
	                type 		= str, 
	                nargs 		= '*', 
					default 	= ['Mercury', 'Venus', 'Earth barycenter', 
					              'Mars barycenter'],
					help		= 'List of bodies for the model, that are '		\
								  'added to the mass of the sun. Make sure to '	\
								  'use planets plus their satellite systems '	\
								  'i.e. "Earth barycenter"')
parser.add_argument('-u', '--units', 
	                action		= 'store', 
	                type 		= str, 
	                nargs 		= '*', 
					default 	= ['au', 'yr', 'msun'],
					help		= 'List of units used in the model. Default is '\
								  'au, yr, and msun')			
parser.add_argument('-s', '--start_time', 
					action		= 'store', 
					type		= str, 
					default 	= '2018-03-30 12:00', 
					help		= 'Date and time of ephemeris in the format '	\
								  '"YYYY-MM-DD HH:MM" the default is set '	\
								  'to 2018-03-30 12:00:00.0 "-s now" sets '		\
								  'start_time to now. If the date entered is '	\
								  'not valid or not in the right format the '	\
								  'date is set to the default value.')
parser.add_argument('-d', '--dt', 
					action		= 'store', 
					type		= float, 
					default 	= 1/365.25/24,	#1 hour 
					help		= 'Timestep for model. Depending on the '		\
								  'underlying integrator dt should be between '	\
								  '10% and 25% of the orbital period or about '	\
								  '0.6rad to 1.5rad. the default here is 1hour'	\
								  'whicgh needs to be set accordingly. ias15 '	\
								  'is using an adaptive timestep.')
parser.add_argument('-t', '--tmax', 
					action		= 'store', 
					type		= float, 
					default 	= 1e4,			#10,000 year
					help		= 'Max duration of model. Default is 10,000yr.')
parser.add_argument('-I', '--integrator', 
					action		= 'store', 
					type		= str, 
					default 	= 'ias15', 
					help		= 'Set integrator used in rebound. Choices are '\
							      '"ias15" (def.), "whfast", "sei", "leapfrog",'\
							      '"hermes", "janus", "mercurius", "bs". See '	\
							      'rebound docs for details.')
parser.add_argument('-G', '--gravity', 
					action		= 'store', 
					type		= str, 
					default 	= 'basic', 
					help		= 'Set gavity module used in rebound. Choices '	\
							      'are "none", "basic" (def.), "compensated" '	\
							      '"tree". See rebound docs for details')
parser.add_argument('-C', '--collision', 
					action		='store',
					type		= str, 
					default 	= 'none',
					help		= 'Set collision module used in rebound. Choices '\
							      'are "none" (default), "direct", "tree", '	\
							      '"mercurius", "direct". See rebound docs for '\
							      'details')
parser.add_argument('-V', '--version', 
					action		= 'version', 
					version		= '%(prog)s 0.1')

clarg = parser.parse_args()

# check and clean up start time 
if clarg.start_time=='now':clarg.start_time = f"{datetime.datetime.now():%Y-%m-%d %H:%M}"
if not addons.validateDate(clarg.start_time):
	print("Runtime Warning: Invalid date format. Date set to 2018-03-30 12:00")
	clarg.start_time = '2018-03-30 12:00'
time = Time(clarg.start_time, scale='utc')

if not clarg.rebuild:

	# check if configfile exists and open it
	CONFIGFILE = clarg.project_name+'.cfg'
	
	if Path(CONFIGFILE).exists():
		cf = configparser.ConfigParser()
		cf.read(CONFIGFILE)
	else:
		clarg.new = True
	
	# check if -n flag is set and write a complete new file now based on commandline
	# arguments and defaults for variables where no command line argument is set 
	# If -n is not set update configfile with new commandline arguments
	if not (clarg.new):
		if (len(sys.argv)) > 2:
			if '-l' in sys.argv:
				cf.set('Bodies', 'list_of_bodies', clarg.list_of_bodies)
			if '-a' in sys.argv:
				cf.set('Bodies', 'add_to_sun', clarg.add_to_sun)
			if '-u' in sys.argv:
				cf.set('Simulation', 'units', clarg.units)
			if '-s' in sys.argv:
				cf.set('Simulation', 'start_time', clarg.start_time) 	 
			if '-s' in sys.argv:
				cf.set('Simulation', 'start_time_jd', str(time.jd))
			if '-d' in sys.argv:
				cf.set('Simulation', 'dt', clarg.dt)
			if '-t' in sys.argv:
				cf.set('Simulation', 'tmax', clarg.tmax)
			if '-I' in sys.argv:
				cf.set('Integrator', 'integrator', clarg.integrator)
			if '-G' in sys.argv:
				cf.set('Integrator', 'gravity', clarg.gravity)
			if '-C' in sys.argv:
				cf.st('Integrator', 'collision', clarg.collision)
			with open(CONFIGFILE, 'w') as configfile:
				cf.write(configfile)
			clarg.new = True	
	else:addons.configWrite(CONFIGFILE, clarg)
else:clarg.new = True

print(clarg.rebuild, clarg.new)
	
BINFILE = clarg.project_name+'-IC.bin'

# Check if configuration file was changed and build a new binfile if it has 
if clarg.new:addons.initReboundBinaryFile(clarg.project_name)

# load simulation from binary file
sim = rb.Simulation.from_file(BINFILE)

sim.status()