import os
import subprocess
import itertools

def replace_block_in_file(filepath, new_content):
    start_marker = "# PYTHON REPLACE START"
    stop_marker = "# PYTHON REPLACE STOP"

    new_lines = []
    inside_block = False
    found_start = False
    found_stop = False

    with open(filepath, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if start_marker in line:
            new_lines.append(line)  # Keep the marker
            new_lines.extend([n + '\n' for n in new_content]) # Add new content
            inside_block = True
            found_start = True
        elif stop_marker in line:
            new_lines.append(line)  # Keep the marker
            inside_block = False
            found_stop = True
        elif not inside_block:
            new_lines.append(line)  # Keep lines outside the markers

    if not (found_start and found_stop):
        print("Error: Could not find both markers in the file.")
        return

    with open(filepath, 'w') as f:
        f.writelines(new_lines)


os.system("cp ../../src/integrator_whfast512.s .")

combinations = []

for stages in range(1,4):
    combinations += list(itertools.product(["newton", "halley"], repeat=stages))
orders = [19, 15, 11, 7]

combinations2 = []
for c in combinations:

    orders2 = list(itertools.product(orders, repeat=len(c)+1))
    for o in orders2:
        c2 = []
        for i,s in enumerate(c):
            if s =="newton":
                c2.append("mm_stiefel_Gs13_avx512 %d"%o[i])
            elif s=="halley":
                c2.append("mm_stiefel_Gs03_avx512 %d"%o[i])
            c2.append(s)

        c2.append("mm_stiefel_Gs13_avx512 %d"%o[-1])
        combinations2.append(c2)



with open("out.txt", "a") as f:
    for i, solvers in enumerate(combinations2):
        print(i,"/",len(combinations2), solvers)
        replace_block_in_file("integrator_whfast512.s", solvers)

        os.system("as -g -o integrator_whfast512.asm_o integrator_whfast512.s")
        os.system("cc -march=native -O3  -std=c99 -Wpointer-arith -D_GNU_SOURCE -fPIC -Wall -g  -Wno-unknown-pragmas -std=c99 -Wpointer-arith -D_GNU_SOURCE -fPIC -Wall -g  -Wno-unknown-pragmas -shared ../../src/*.o -lm -lrt integrator_whfast512.asm_o -o librebound.so")
        output = subprocess.check_output(["./rebound"], text=True)

        print(",".join(solvers) + "\t"+ output.strip(), file=f, flush=True)
        print(output.strip(), flush=True)
