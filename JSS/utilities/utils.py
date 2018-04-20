import re
import rebound as rb
import datetime
import math

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


def initReboundBinaryFile():
###############################################################################
# This subroutine is still very static, designed to cater for JSS-LR-8S and
# needs to be changed in order to be more flexible (i.e. Configiuration file)
# containing:
#   List of bodies
#   List of bodies to add to the Sun
#   Model units
#   time of ephemeris (gregorian or julian. astropy?)
#   Name of project
# Should not be too hard to implement ....
###############################################################################

# define variables pretending to be constants

    MASS_SOLAR          = getMass('Sun')              # in kg
    MASS_MERCURY        = getMass('Mercury')          # in kg
    MASS_VENUS          = getMass('Venus')            # in kg
    MASS_EARTH          = getMass('Earth barycenter') # in kg (Earth+Moon)
    MASS_MARS           = getMass('Mars barycenter')  # in kg (Mars+Moons)
    
# Set bodies for model

    bodies  =  ["Jupiter", "Metis", "Adrastea", "Amalthea", "Thebe", "Io", "Europa", 
                "Ganymede", "Callisto", "Sun", "Saturn barycenter", "Uranus barycenter", 
                "Neptune barycenter"]
    
    for i in range(0,len(bodies),1):bodies[i] = str(getNAIF(bodies[i]))
    
# Set rebound integrator conditions

    model               = rb.Simulation()
    model.units         = ("yr", "au", "Msun")  # chosen yr instead of yr/2pi
    model.integrator    = "whfast"              # whfast integrator (no collision)
    #model.integrator   = "IAS15"               # IAS15 integrator (collision)
    time_of_ephemeris   = 2458208.              # Julian Day of 30/03/2018 12:00
    model.t             = time_of_ephemeris
    bin_file            = '../data/JSS-LR-8S-InCond-JD'+str(time_of_ephemeris)+'.bin'
    
    model.add(bodies, date='2018-03-30 12:00')  # Julian day 2458208.
    
# Set particle hash to address particles by their name rather than their index

    for i in range(0,len(model.particles),1):
        model.particles[i].hash = getNAIF(int(bodies[i]))
    
# Throw the inner planets into the Sun "Boohahaha!"

    model.particles['Sun'].m        = (MASS_SOLAR+      \
                                       MASS_MERCURY+    \
                                       MASS_VENUS+      \
                                       MASS_EARTH+      \
                                       MASS_MARS)/      \
                                       MASS_SOLAR
    
# Save initial conditions to binary file

    model.status()
    model.save(bin_file)
    