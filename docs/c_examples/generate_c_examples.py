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
                fd.write("\n\nThis example is located in the directory `examples/"+problemc.split("/")[2]+"`\n\n")
                if did_output==0:
                    print("Warning: Did not find description in "+problemc)
    print("Converted %d C examples."%count)
if __name__ == "__main__":
    run()
