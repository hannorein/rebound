import configparser

cfin = configparser.ConfigParser()

cfin.read('example.ini')
list_of_bodies 			= cfin['Bodies']['list_of_bodies']
add_to_sun				= cfin['Bodies']['add_to_sun']



units					= cfin['Simulation']['units']


list_of_bodies = configStr2List(list_of_bodies)
add_to_sun = configStr2List(add_to_sun)
units =	configStr2List(units)

print(list_of_bodies, list_of_bodies[1])
print(add_to_sun, add_to_sun[1])
print(units, units[1])

def configStr2List(convert_me):
	remove = "'[]"
	for char in remove:convert_me = convert_me.replace(char,'')
	char = ", "
	convert_me = convert_me.replace(char,',')
	convert_me = convert_me.split(',')
	for i in range(len(convert_me)):convert_me[i].lstrip