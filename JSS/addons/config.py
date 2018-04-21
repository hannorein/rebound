import addons
from argparse import ArgumentParser as ap
import configparser
from astropy.time import Time
import datetime
import sys


parser = ap(description='Reads command line arguments and creates configuration'\
			' file. Start with ctrebound <project_name> where <project_name> '	\
			'is required. With no other arguments provided, the configuration '	\
			'file for <project_name> is used. If there is no configuration '	\
			'file in the data directory, the program will create a new config '	\
			'file using the defaults as well as arguments if provided. -o will '\
		    'overwrite an existing file without warning and -k will update an '	\
		    'existing configuration file only with the arguments provided' )
parser.add_argument('project_name', 
					type		= str,
					help		= 'Id of the project. (example: JSS-LR-12S for '\
								  '"Jupiter-Satellite-System-Long-Run-12-Satell'\
								  'ites). <project_name> is used in all project'\
								  ' related files.  Use whatever rocks your '	\
								  'boat! ')
parser.add_argument('-k', '--keep', 											
					action		= 'store_const',
					const 		= True,
					help 		= 'Keep and update configuartion file with the '\
								  'arguments provided in the command line.')			
parser.add_argument('-o', '--overwrite', 											
					action		= 'store_const',
					const 		= True,
					help 		= 'Overwrite the configuartion file with the '	\
								  'arguments provided in the command line.')
parser.add_argument('-l', '--list_of_bodies', 
	                action		= 'store', 
	                type 		= str, 
	                nargs 		= '*', 
					default 	= ['Jupiter', 'Metis', 'Adrestea', 'Amalthea',
					              'Thebe', 'Io', 'Europa', 'Ganymede', 'Calisto',
					              'Sun', 'Saturn barycenter', 'Uranus barycenter', 
					              'Neptune barycenter'],
					help		= 'List of bodies for the model, separeded by '	\
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
					default 	= '2018-03-30 12:00:00.0', 
					help		= 'Date and time of ephemeris in the format '	\
								  '"YYYY-MM-DD HH:MM:00.0" the default is set '	\
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
parser.add_argument('-B', '--boundary', 
					action		='store',
					type		= str, 
					default 	= 'none',
					help		= 'Set boundary module used in rebound. Choices '\
							      'are "none" (default), "open", "periodic", '	\
							      '"shear". See rebound docs for details')
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

print(sys.argv[1:])

if clarg.start_time=='now':clarg.start_time = str(datetime.datetime.now())
if not addons.validateDate(clarg.start_time):
	clarg.start_time 			= '2018-03-30 12:00:00.0'

time 							= Time(clarg.start_time, scale='utc')

cfout = configparser.ConfigParser()
cfout['Bodies'] 				= {'list_of_bodies'	: clarg.list_of_bodies,
					   			   'add_to_sun'		: clarg.add_to_sun}
cfout['Simulation'] 			= {'project_name'	: clarg.project_name,
								   'units'			: clarg.units,
								   'start_time'		: clarg.start_time,
								   'start_time_jd'	: time.jd,
								   'dt'				: clarg.dt,
								   'tmax'			: clarg.tmax}
cfout['Integrator'] 			= {'integrator'		: clarg.integrator,
								   'gravity'		: clarg.gravity,
								   'boundary'		: clarg.boundary,
								   'collision'		: clarg.collision}

with open(clarg.project_name+'.ini', 'w') as configfile:
	cfout.write(configfile)