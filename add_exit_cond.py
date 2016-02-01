#!/usr/local/bin/python

import sys

if len(sys.argv) is not 6:
    print('Error\n*****\nMust pass (in order) name of parameter, name of type in c (e.g. "struct reb_particle*"), the ctype type (e.g., "POINTER(Particle)"), the number code for the enum (shouldnt conflict with codes at the top of src/rebound.h in REB_STATUS, e.g., 736), and a description to be shown in python when error occurs (e.g. "a particle got too close to particle blah")')
    sys.exit()

param_name = str(sys.argv[1])
param_type_c = str(sys.argv[2])
param_ctype = str(sys.argv[3])
param_enum_code = str(sys.argv[4])
description = str(sys.argv[5])

exception = param_name.capitalize()

tab = "    " # 4 space tab

with open("src/rebound.h") as f:
    rebh = f.readlines()

with open("src/rebound.h", "w") as f:
    for i in range(len(rebh)):
        f.write(rebh[i])
        if "struct reb_simulation {" in rebh[i]:
            f.write(tab+"{0} {1};\n".format(param_type_c, param_name))
        if "enum REB_STATUS {" in rebh[i]:
            f.write(tab+"REB_{0} = {1},\n".format(param_name.upper(), param_enum_code))

with open("src/rebound.c") as f:
    rebc = f.readlines()

with open("src/rebound.c", "w") as f:
    for i in range(len(rebc)):
        f.write(rebc[i])
        if "if (r->heartbeat){ r->heartbeat(r); }" in rebc[i]:
            f.write(tab+"if (r->{0}){{\n".format(param_name))
            f.write(tab+tab+"// Check for custom condition\n")
            f.write(tab+tab+"if (fill_in_here){\n")
            f.write(tab+tab+tab+"r->status = REB_{0};\n".format(param_name.upper()))
            f.write(tab+tab+"}\n")
            f.write(tab+"}\n")
        if "void reb_init_simulation(" in rebc[i]:
            f.write("    r->{0} = 0;\n".format(param_name))

with open("rebound/simulation.py") as f:
    sim = f.readlines()

with open("rebound/simulation.py", "w") as f:
    for i in range(len(sim)):
        if "from . import clibrebound" in sim[i]:
            f.write(sim[i].rstrip()+", {0}\n".format(exception))
        else:
            f.write(sim[i])
        if "Simulation._fields_ = " in sim[i]:
            f.write(tab+tab+tab+tab+'("{0}", {1}),\n'.format(param_name, param_ctype))
        if "ret_value = clibrebound.reb_integrate" in sim[i]:
            f.write(tab+tab+tab+"if ret_value == {0}:\n".format(param_enum_code))
            f.write(tab+tab+tab+tab+'raise {0}("{1}")\n'.format(exception, description))

with open("rebound/__init__.py") as f:
    init = f.readlines()

with open("rebound/__init__.py", "w") as f:
    for i in range(len(init)):
        if "__all__ = " in init[i]:
            init[i] = init[i].rstrip()
            f.write(init[i][:-1] + ', "{0}"]\n'.format(exception))
        else:
            f.write(init[i])
        if "# Exceptions" in init[i]:
            f.write("class {0}(Exception):\n".format(exception))
            f.write(tab+"pass\n")
