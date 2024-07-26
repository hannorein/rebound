#!python
# This script automatically creates a list of examples by reading the header in all problem.c files.
import subprocess
ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()

with open("version.txt") as f:
    reboundversion = f.readlines()[0].strip()
    print("Updating version to "+reboundversion)

with open("README.md") as f:
    readme = f.readlines()

with open("README.md","w") as f:
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

with open("web_client/shell_rebound_webgl.html") as f:
    reboundlines = f.readlines()
    for i,l in enumerate(reboundlines):
        if "<!-- VERSIONLINE -->" in l:
            reboundlines[i] = "                  REBOUND v" + reboundversion + "  <!-- VERSIONLINE -->\n"

    with open("web_client/shell_rebound_webgl.html", "w") as f:
        f.writelines(reboundlines)

with open("web_client/shell_rebound_console.html") as f:
    reboundlines = f.readlines()
    for i,l in enumerate(reboundlines):
        if "<!-- VERSIONLINE -->" in l:
            reboundlines[i] = "                  REBOUND v" + reboundversion + "  <!-- VERSIONLINE -->\n"

    with open("web_client/shell_rebound_console.html", "w") as f:
        f.writelines(reboundlines)

with open("web_client/shell_rebound.html") as f:
    reboundlines = f.readlines()
    for i,l in enumerate(reboundlines):
        if "<!-- VERSIONLINE -->" in l:
            reboundlines[i] = "                  REBOUND v" + reboundversion + "  <!-- VERSIONLINE -->\n"

    with open("web_client/shell_rebound.html", "w") as f:
        f.writelines(reboundlines)

shortversion = reboundversion
while shortversion[-1] != '.':
    shortversion = shortversion[:-1]
    
shortversion = shortversion[:-1]

# find changelog
with open("changelog.md") as f:
    found_start = 0
    changelog = ""
    cl = f.readlines()
    for l in cl:
        if found_start == 0 and l.startswith("### Version"):
            if reboundversion in l:
                found_start = 1
                continue
        if found_start == 1 and l.startswith("### Version"):
            found_start = 2
        if found_start == 1:
            changelog += l

if found_start != 2 or len(changelog.strip())<5:
    raise RuntimeError("Changelog not found")

with open("_changelog.tmp", "w") as f:
    f.writelines(changelog.strip()+"\n")

print("----")
print("Changelog:\n")
print(changelog.strip())
print("----")
print("Next:")
print("\ngit commit -a -m \"Updating version to "+reboundversion+"\"")
print("git tag "+reboundversion+" && git push && git push --tags")
print("gh release create "+reboundversion+" --notes-file _changelog.tmp")
print("----")
print("Might also want to push a rebound.html to this release:")
print("gh release upload "+reboundversion+" web_client/rebound.html")
