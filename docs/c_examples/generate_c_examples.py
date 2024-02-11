# -*- coding: utf-8 -*-
import glob

# C Example update
def run(*args, **kwargs):
    print("Generating C examples.")
    count = 0
    for problemc in glob.glob("examples/*/problem.c"):
        count += 1
        cname = problemc.split("/")[1]
        with open("docs/c_examples/"+cname+".md","w") as fd:
            will_output = 0
            livepreview=1
            # Manual exception for file viewer
            if "_viewer" in cname:
                livepreview=0
            if "screenshots" in cname:
                livepreview=0
            try:
                with open("examples/"+cname+"/Makefile","r") as mfd:
                    Makefile = mfd.read()
                    if "export MPI=1" in Makefile:
                        livepreview=0
                    if "export OPENMP=1" in Makefile:
                        livepreview=0
            except:
                print("Warning: Makefile error in "+problemc)
            
            with open(problemc) as pf:
                did_output=0
                empty_lines = 0
                for line in pf:
                    if line[0:3] == "/**":
                        will_output += 1
                    if line[0:3] == " */":
                        will_output = -1
                        line = ""
                        fd.write("\n\n```c\n")
                    if will_output>1:
                        if will_output == 2:
                            line = "   # "+line[3:].strip() + " (C)\n"
                            if livepreview == 1:
                                line += "!!! example \"Try it out this example!\"\n"
                                line += "    REBOUND has been compiled with emscripten to WebAssembly.\n"
                                line += "    This lets you run this example interactively from within your browser at almost native speed.\n"
                                line += "    No installation is required.\n"
                                line += "    [Click here](../../emscripten_c_examples/"+cname+"/) to try it out.\n"
                        will_output = 2
                        if len(line[3:].strip())==0:
                            fd.write("\n\n"+line[3:].strip())
                        else:
                            fd.write(line[3:].strip() + " " )
                    if will_output==-1:
                        fd.write("" +line.rstrip() + "\n" )
                        did_output = 1
                    if will_output>0:
                        will_output += 1
                fd.write("```\n")
                fd.write("\n\nThis example is located in the directory `examples/"+problemc.split("/")[1]+"`\n\n")
                if did_output==0:
                    print("Warning: Did not find description in "+problemc)
    print("Converted %d C examples."%count)
if __name__ == "__main__":
    run()
