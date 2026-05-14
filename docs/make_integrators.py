import rebound
import ctypes
import os
from rebound.binarydata_field_descriptor import REB_BINARYDATA_DTYPE
os.makedirs("integrators", exist_ok=True)

# Get all integrators
clibrebound = rebound.clibrebound
clibrebound.reb_integrators_registered.restype = ctypes.POINTER(ctypes.c_char_p)
res = clibrebound.reb_integrators_registered()
i = 0
names = []
while True:
    if res[i] == None:
        break
    else:
        names.append(res[i].decode("ascii"))
    i = i+1
clibrebound.reb_free(res)

for name in names:
    print(name)
    with open("integrators/"+name+".md", "w") as file:
        sim = rebound.Simulation()
        sim.integrator = name
        file.write("# " + name.upper()+"\n\n")
        file.write(sim.integrator.callbacks.documentation.decode("utf-8") + "\n\n")
        fdlist = sim.integrator.callbacks.field_descriptor_list
        first_attr = True
        for fd in fdlist:
            fd_doc = fd.doc()
            if fd_doc:
                if first_attr:
                    file.write("## Attributes \n\n")
                    first_attr = False

                if fd.dtype not in REB_BINARYDATA_DTYPE:
                    dtype = "unknown datatype"
                else:
                    dtype = str(REB_BINARYDATA_DTYPE[fd.dtype].__name__)
                    if dtype.startswith("c_"):
                        dtype = dtype[2:]
                    if dtype.endswith("_p"):
                        dtype = dtype[:-2] + "*"
                    dtype = "`"+dtype+"`"
                file.write("### " + fd.name.decode("utf-8") + " (" +dtype+ ") \n")
                file.write(fd.documentation.decode("utf-8")+"\n\n")
                edl = fd.enum_descriptor_list
                if edl:
                    file.write("#### Supported values\n\n")
                    file.write("| enum constant | identifier |\n")
                    file.write("| ------------- | ---------- |\n")
                    for ed in edl:
                        file.write("| "+ str(ed.value) + " | `"+ed.name.decode("utf-8") + "` |\n")
                file.write("\n")



with open("."+"./mkdocs.yml", "r") as file:
    lines = file.readlines()
    
output_lines = []
inside_markers = False
replaced = False

# 2. Process the file line by line
for line in lines:
    if line.strip() == "# Integrators --START--":
        inside_markers = True
        output_lines.append(line)
        continue
        
    if line.strip() == "# Integrators --END--":
        inside_markers = False
        for name in names:
            output_lines.append("      - integrators/"+name+".md\n")
        output_lines.append(line)
        replaced = True
        continue
    if not inside_markers:
        output_lines.append(line)

with open("."+"./mkdocs.yml", "w") as file:
    file.write("".join(output_lines))


