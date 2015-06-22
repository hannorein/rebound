import json
import sys
print "import matplotlib as mpl"
print "mpl.use(\"Agg\")"

if len(sys.argv)!=2:
    print "Usage: ipynb2py.py FILENAME"
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
print code
