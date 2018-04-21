import configparser
import addons

cfin = configparser.ConfigParser()

cfin.read('example.ini')
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