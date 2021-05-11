#!/usr/bin/python
# This script automatically creates a list of examples by reading the header in all problem.c files.
import glob
import subprocess
ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()

with open("version.txt") as f:
    reboundversion = f.readlines()[0].strip()
    print("Updating version to "+reboundversion)

with open("README.md") as f:
    readme = f.readlines()

keep_lines_after_header = 5

with open("README.md","w") as f:
    start_delete = -1
    for i in range(0,len(readme)):
        # [![Version](https://img.shields.io/badge/rebound-v3.17.0-green.svg?style=flat)](https://rebound.readthedocs.org)
        if "![Version]" in readme[i]:
            readme[i] = "[![Version](https://img.shields.io/badge/rebound-v"+reboundversion+"-green.svg?style=flat)](https://rebound.readthedocs.org)\n"
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
        if "GITHASHAUTOUPDATE" in l:
            setuplines[i] = "    ghash_arg = \"-DGITHASH="+ghash+"\" #GITHASHAUTOUPDATE\n"

    with open("setup.py", "w") as f:
        f.writelines(setuplines)

shortversion = reboundversion
while shortversion[-1] != '.':
    shortversion = shortversion[:-1]
    
shortversion = shortversion[:-1]

print("To commit, copy and paste:")
print("\ngit commit -a -m \"Updating version to "+reboundversion+"\"")
