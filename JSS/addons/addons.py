import re
import rebound as rb
import datetime
import math
import configparser
from astropy.time import Time

def configWrite(file_name, clarg, sys_argv, update=True):
    # clarg = command line arguments from argparser
    if update:
        cf = configparser.ConfigParser()
        cf.read(file_name)
    else:
        cf = configparser.ConfigParser()
        cf.read('.rebounddefault.cfg')
    
    if '-l' in sys_argv:
        cf.set('Bodies', 'list_of_bodies', str(clarg.list_of_bodies))
    if '-a' in sys_argv:
        cf.set('Bodies', 'add_to_sun', str(clarg.add_to_sun))
    if '-u' in sys_argv:
        cf.set('Simulation', 'units', str(clarg.units))
    if '-s' in sys_argv:
        if clarg.start_time=='now':
            clarg.start_time = f"{datetime.datetime.now():%Y-%m-%d %H:%M}"
        if not validateDate(clarg.start_time):
            print("Warning: Invalid date format. Date set to 2000-01-01 12:00")
            clarg.start_time = '2000-01-01 12:00'
        time = Time(clarg.start_time, scale='utc')
        cf.set('Simulation', 'start_time', clarg.start_time)
        cf.set('Simulation', 'start_time_jd', str(time.jd))
    if '-d' in sys_argv:
        cf.set('Simulation', 'dt', clarg.dt)
    if '-t' in sys_argv:
        cf.set('Simulation', 'tmax', clarg.tmax)
    if '-I' in sys_argv:
        cf.set('Integrator', 'integrator', clarg.integrator)
    if '-G' in sys_argv:
        cf.set('Integrator', 'gravity', clarg.gravity)
    if '-C' in sys_argv:
        cf.set('Integrator', 'collision', clarg.collision)
    with open(file_name, 'w') as configfile:
        cf.write(configfile)

    return 0

def validateDate(date_text):
    # checks if date is in format YYYY-MM-DD HH:MM:SS.MS
    try:
        datetime.datetime.strptime(date_text, '%Y-%m-%d %H:%M')
        return True
    except ValueError:
        return False

def configStr2List(convert_me):
    remove = "'[]"
    for char in remove:convert_me = convert_me.replace(char,'')
    char = ", "
    convert_me = convert_me.replace(char,',')
    convert_me = convert_me.split(',')
    for i in range(len(convert_me)):convert_me[i].lstrip
    return convert_me

def getMass(body=None):
	if type(body) is str:
		body = rb.horizons.getNAIF(body)
	elif type(body) is int:
		body = str(body)
	else:
		raise AttributeError("Body needs to be a name or NAIF identifier.")

	mass = float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int(body), rb.horizons.HORIZONS_MASS_DATA).group(1))
	mass /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)
	return mass

def initReboundBinaryFile(project_name, sys_argv):

    BINFILE = project_name+'-InitialConditions.bin'
    CONFIGFILE = project_name+'.cfg'
    REBUILD_FLAGS = ['-n', '-r', '-l', '-a', '-u', '-s']

# read in config file

    cfin = configparser.ConfigParser()
    cfin.read(CONFIGFILE)

# Set up particles for model and particles of which the mass will be added to the Sun
    particles       = configStr2List(cfin['Bodies']['list_of_bodies'])
    add_to_sun      = configStr2List(cfin['Bodies']['add_to_sun'])
    if add_to_sun[0]:particles = particles + add_to_sun

# Set rebound integrator conditions
    if [i for i in REBUILD_FLAGS if i in sys_argv]:
        sim         = rb.Simulation()
        sim.units   = configStr2List(cfin['Simulation']['units'])

        sim.add(particles, date=cfin['Simulation']['start_time'])

        if add_to_sun[0]:
            for particle in particles:
                if particle in add_to_sun: 
                    sim.particles['Sun'].m += sim.particles[particle].m
                    sim.remove(hash=particle)
    else:
        sim         = rb.Simulation.from_file(BINFILE)

    sim.integrator  = cfin['Integrator']['integrator']
    sim.t           = float(cfin['Simulation']['start_time_jd'])
    sim.dt          = float(cfin['Simulation']['dt'])
    sim.gravity     = cfin['Integrator']['gravity']
    sim.collision   = cfin['Integrator']['collision']
    tmax            = float(cfin['Simulation']['tmax'])


# Save initial conditions to binary file
    sim.save(BINFILE)

    return 0
    