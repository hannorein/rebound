import json
import sys
exec("import matplotlib as mpl")
exec("mpl.use(\"Agg\")")

if len(sys.argv)!=2:
    print("Usage: ipynb2py.py FILENAME")
    exit(1)
with open(sys.argv[1]) as data_file:    
    ipynb = json.load(data_file)

code = ""
for c in ipynb["cells"]:
    if c["cell_type"] == "code":
        source = c["source"]
        for s in source:
            if s[0] != "%":
                code += s.rstrip('\n')+"\n"
import socket
try:
    exec(code)
except socket.error: 
    print("A socket error occured. This is most likely due to a timeout in the NASA Horizons connections. We catch this exception here and ignore is.")
    pass
