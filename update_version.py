#!/usr/bin/python
# This script automatically creates a list of examples by reading the header in all problem.c files.
import glob

with open("version.txt") as f:
    reboundversion = f.readlines()[0].strip()
    print("Updating version to "+reboundversion)

with open("README.rst") as f:
    readme = f.readlines()

keep_lines_after_header = 5

with open("README.rst","w") as f:
    start_delete = -1
    for i in xrange(0,len(readme)):
        if "badge/rebound-v" in readme[i]:
            readme[i] = ".. image:: http://img.shields.io/badge/rebound-v"+reboundversion+"-green.svg?style=flat\n"
        f.write(readme[i])

with open("src/rebound.c") as f:
    reboundlines = f.readlines()
    for i,l in enumerate(reboundlines):
        if "**VERSIONLINE**" in l:
            reboundlines[i] = "const char* reb_version_str = \""+reboundversion+"\";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.\n"

    with open("src/rebound.c", "w") as f:
        f.writelines(reboundlines)

with open("setup.py") as f:
    setuplines = f.readlines()
    for i,l in enumerate(setuplines):
        if "version='" in l:
            setuplines[i] = "    version='"+reboundversion+"',\n"

    with open("setup.py", "w") as f:
        f.writelines(setuplines)

shortversion = reboundversion
while shortversion[-1] is not '.':
    shortversion = shortversion[:-1]
    
shortversion = shortversion[:-1]

with open("doc/conf.py") as f:
    conflines = f.readlines()
    for i,l  in enumerate(conflines):
        if "version =" in l:
            conflines[i] = "version = '"+shortversion+"'\n"
        if "release =" in l:
            conflines[i] = "release = '"+reboundversion+"'\n"

    with open("doc/conf.py", "w") as f:
        f.writelines(conflines)
print("To commit, copy and paste:")
print("\ngit commit -a -m \"Updating version to "+reboundversion+"\"")
