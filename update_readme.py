#!/usr/bin/python
# This script automatically creates a list of examples by reading the header in all problem.c files.
import glob

with open("README.md") as f:
    readme = f.readlines()

keep_lines_after_header = 5
with open("README.md","w") as f:
    start_delete = -1
    for i in xrange(0,len(readme)):
        if readme[i]=="Examples\n":
            if readme[i+1]=="-----------------------\n":
                start_delete = i + keep_lines_after_header
        if start_delete!=-1 and readme[i+2]=="-----------------------\n":
            for problemc in glob.glob("./examples/*/problem.c"):
                will_output = 0
                with open(problemc) as pf:
                    did_output=0
                    for line in pf:
                        if line[0:10] == " * @detail":
                            will_output = 1
                            did_output = 1
                            line = " * " + line[10:]
                            f.write("\n*  **examples/"+problemc.split("/")[2]+"**\n\n")
                            f.write("  This example is using the following modules:  \n");
                            with open("examples/"+problemc.split("/")[2]+"/Makefile") as mf:
                                for linem in mf:
                                    if linem.strip()[0:6] == "ln -fs":
                                        sourcefile = linem.strip().split(" ")[2].split("/")[-1]
                                        if sourcefile!="problem.c":
                                            f.write("  `"+sourcefile+"`\n")
                            f.write("\n");
                    
                        if line[0:11] == " * @section":
                            will_output = 0
                        if will_output==1:
                            f.write("  "+line[3:].strip()+"\n")
                    if did_output==0:
                        print "Warning: Did not find description in "+problemc
            start_delete = -1
        if start_delete==-1 or i<start_delete:
            f.write(readme[i])

