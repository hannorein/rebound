import configparser
import addons

cfin = configparser.ConfigParser()

cfout.read('test.ini')
cfout.set('Bodies', 'list_of_bodies', 'Ha I changed you again!')
with open('test.ini', 'w') as configfile:
	cfout.write(configfile)



cfin.read('JSS-LR-8S.ini')
list_of_bodies 			= addons.configStr2List(cfin['Bodies']['list_of_bodies'])
add_to_sun				= addons.configStr2List(cfin['Bodies']['add_to_sun'])
project_name			= cfin['Simulation']['project_name']
units					= addons.configStr2List(cfin['Simulation']['units'])
start_time 				= cfin['Simulation']['start_time']
start_time_jd 			= cfin['Simulation']['start_time_jd']
dt 						= cfin['Simulation']['dt']
tmax 					= cfin['Simulation']['tmax']
integrator 				= cfin['Integrator']['integrator']
gravity  				= cfin['Integrator']['gravity']
boundary  				= cfin['Integrator']['boundary']
collision  				= cfin['Integrator']['collision']


cfout = configparser.ConfigParser()
cfout['Bodies'] 				= {'list_of_bodies'	: 'clarg.list_of_bodies',
					   			   'add_to_sun'		: 'clarg.add_to_sun'}
cfout['Simulation'] 			= {'project_name'	: 'clarg.project_name',
								   'units'			: 'clarg.units',
								   'start_time'		: 'clarg.start_time',
								   'start_time_jd'	: 'time.jd',
								   'dt'				: 'clarg.dt',
								   'tmax'			: 'clarg.tmax'}
cfout['Integrator'] 			= {'integrator'		: 'clarg.integrator',
								   'gravity'		: 'clarg.gravity',
								   'boundary'		: 'clarg.boundary',
								   'collision'		: 'clarg.collision'}
with open('test.ini', 'w') as configfile:
	cfout.write(configfile)