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
        cf.set('Bodies', 'list_of_bodies', clarg.list_of_bodies)
    if '-a' in sys_argv:
        cf.set('Bodies', 'add_to_sun', clarg.add_to_sun)
    if '-u' in sys_argv:
        cf.set('Simulation', 'units', clarg.units)
    if '-s' in sys_argv:
        if clarg.start_time=='now':
            clarg.start_time = f"{datetime.datetime.now():%Y-%m-%d %H:%M}"
        if not validateDate(clarg.start_time):
            print("Runtime Warning: Invalid date format. Date set to 2018-03-30 12:00")
            clarg.start_time = '2018-03-30 12:00'
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
        cf.st('Integrator', 'collision', clarg.collision)
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
		body = getNAIF(body)
	elif type(body) is int:
		body = str(body)
	else:
		raise AttributeError("Body needs to be a name or NAIF identifier.")

	mass = float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int(body), rb.horizons.HORIZONS_MASS_DATA).group(1))
	mass /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)
	return mass

def getNAIF(name_or_id=None):
    
    NAIF_NAME_OR_ID={
        "Solar system barycenter":  0,
        "Mercury barycenter":       1,
        "Venus barycenter":         2,
        "Earth barycenter":         3,
        "Mars barycenter":          4,
        "Jupiter barycenter":       5,
        "Saturn barycenter":        6,
        "Uranus barycenter":        7,
        "Neptune barycenter":       8,
        "Pluto barycenter":         9,
        "Sun":                      10,
        "Mercury":                  199,
        "Venus":                    299,
        "Earth":                    399,
        "Mars":                     499,
        "Jupiter":                  599,
        "Saturn":                   699,
        "Uranus":                   799,
        "Neptune":                  899,
        "Pluto":                    999,
        "Moon":                     301,
        "Phobos":                   401,
        "Deimos":                   402,
        "Io":                       501,
        "Europa":                   502,
        "Ganymede":                 503,
        "Callisto":                 504,
        "Amalthea":                 505,
        "Himalia":                  506,
        "Elara":                    507,
        "Pasiphae":                 508,
        "Sinope":                   509,
        "Lysithea":                 510,
        "Carme":                    511,
        "Ananke":                   512,
        "Leda":                     513,
        "Thebe":                    514,
        "Adrastea":                 515,
        "Metis":                    516,
        "Mimas":                    601,
        "Enceladus":                602,
        "Tethys":                   603,
        "Dione":                    604,
        "Rhea":                     605,
        "Titan":                    606,
        "Hyperion":                 607,
        "Iapetus":                  608,
        "Phoebe":                   609,
        "Janus":                    610,
        "Epimetheus":               611,
        "Helene":                   612,
        "Telesto":                  613,
        "Calypso":                  614,
        "Atlas":                    615,
        "Prometheus":               616,
        "Pandora":                  617,
        "Pan":                      618,
        "Ariel":                    701,
        "Umbriel":                  702,
        "Titania":                  703,
        "Oberon":                   704,
        "Miranda":                  705,
        "Cordelia":                 706,
        "Ophelia":                  707,
        "Bianca":                   708,
        "Cressida":                 709,
        "Desdemona":                710,
        "Juliet":                   711,
        "Portia":                   712,
        "Rosalind":                 713,
        "Belinda":                  714,
        "Puck":                     715,
        "Caliban":                  716,
        "Sycorax":                  717,
        "1986U10":                  718,
        "Triton":                   801,
        "Nereid":                   802,
        "Naiad":                    803,
        "Thalassa":                 804,
        "Despina":                  805,
        "Galatea":                  806,
        "Larissa":                  807,
        "Proteus":                  808,
        "Charon":                   901
    }

    if type(name_or_id) is int:
        if name_or_id in NAIF_NAME_OR_ID.values():
            return list(NAIF_NAME_OR_ID.keys())[list(NAIF_NAME_OR_ID.values()).index(name_or_id)]
        else:
            raise AttributeError("Unknown NAIF id. For a list of valid codes and ids call utils.getNAIF('print')")
    elif type(name_or_id) is str:
        if name_or_id in NAIF_NAME_OR_ID:
            return NAIF_NAME_OR_ID[name_or_id]
        elif name_or_id == 'print':
            for key in NAIF_NAME_OR_ID: print("ID: {} = {}".format(NAIF_NAME_OR_ID[key], key))
        else:
            raise AttributeError("Unknown NAIF name. For a list of valid codes and ids call utils.getNAIF('print')")
    else:
        raise AttributeError("NAIF name or id code not listed. For a list of valid codes and ids call util.getNAIF('print')")


def initReboundBinaryFile(project_name):

    BINFILE = project_name+'-IC.bin'

# read in config file
    CONFIGFILE = project_name+'.cfg'

    cfin = configparser.ConfigParser()
    cfin.read(CONFIGFILE)

# Set up bodies for model and get NAIF id numbers for those bodies
    bodies = configStr2List(cfin['Bodies']['list_of_bodies'])
    for i in range(len(bodies)):
        bodies[i] = str(getNAIF(bodies[i]))
    
    add_to_sun = configStr2List(cfin['Bodies']['add_to_sun'])
    for i in range(0,len(add_to_sun),1):
        add_to_sun[i] = str(getNAIF(add_to_sun[i]))

    bodies     = bodies + add_to_sun
    add_to_sun = configStr2List(cfin['Bodies']['add_to_sun'])

# Set rebound integrator conditions
    sim             = rb.Simulation()
    sim.units       = configStr2List(cfin['Simulation']['units'])
    sim.integrator  = cfin['Integrator']['integrator']
    sim.t           = float(cfin['Simulation']['start_time_jd'])
    sim.dt          = float(cfin['Simulation']['dt'])
    sim.gravity     = cfin['Integrator']['gravity']
    sim.collision   = cfin['Integrator']['collision']
    tmax            = float(cfin['Simulation']['tmax'])

    sim.add(bodies, date=cfin['Simulation']['start_time'])

# Set particle hash to address particles by their name rather than their index
# and add up masses of bodies in add_to_sun
    add_mass = 0
    for i in range(0,len(sim.particles),1):
        body = getNAIF(int(bodies[i]))
        sim.particles[i].hash = body
        if body in add_to_sun:add_mass = add_mass+sim.particles[i].m

# Throw the inner planets into the Sun "Boohahaha!"        
    sim.particles['Sun'].m = sim.particles['Sun'].m + add_mass

# Remove bodies in add_to_sun from simulation
    for body in add_to_sun:
        sim.remove(hash=body)

# Save initial conditions to binary file
    sim.save(BINFILE)

    return 0
    