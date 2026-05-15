import re
import os

# C Example update
print("Generating C examples.")
with open("mkdocs.yml","r") as f:
    config = f.readlines()

count = 0
for l in config:                  
    m = re.search(r'.*c_examples/([^\.]+)\.md.*', l)
    if m:
        count += 1
        name = m.group(1)
        print("Working on ", name)
        problemc = "examples/"+name+"/problem.c"                  
        with open("docs/c_examples/"+name+".md","w") as fd:
            will_output = 0
            livepreview=1
            # Manual exception for file viewer
            if "_viewer" in name:
                livepreview=0
            if "screenshots" in name:
                livepreview=0
            try:
                with open("examples/"+name+"/Makefile","r") as mfd:
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
                server_used = 0
                for line in pf:
                    if "reb_simulation_start_server" in line:
                        server_used = 1
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
                                line += "    [Click here](../../emscripten_c_examples/"+name+"/) to try it out.\n"
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
            if livepreview:
                os.makedirs("docs/c_examples_emscription/"+name, exist_ok=True)
                print("Compiling with emscripten..." 
                if server_used == 0:
                    os.system("source ../emsdk/emsdk_env.sh && emcc -O3 -Isrc/ src/*.c $dir/problem.c -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -s -sASYNCIFY -sALLOW_MEMORY_GROWTH -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file web_client/shell_rebound_console.html -o site/emscripten_c_$dir/index.html")
                else:
                    os.system("source ../emsdk/emsdk_env.sh && emcc -O3 -Isrc/ src/*.c $dir/problem.c -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sEXPORTED_RUNTIME_METHODS=\"callMain\" --shell-file web_client/shell_rebound_webgl.html -o site/emscripten_c_$dir/index.html"

print("Converted %d C examples."%count)
