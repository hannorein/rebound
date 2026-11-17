import json
import sys
exec("import matplotlib as mpl")
exec("mpl.use(\"Agg\")")

if len(sys.argv)!=2:
    # Run all ipynb files in this directory except the ones listed below.
    DO_NOT_RUN = ["RealtimeVisualizations.ipynb","Cheartbeat.ipynb","ChaoticHyperion.ipynb", "WHFast512.ipynb"]

    print("Running *all* ipynb files in this directory.")
    import subprocess
    import glob
    for filename in glob.glob("*.ipynb"):
        if filename in DO_NOT_RUN:
            print("Skipping " + filename)
            continue

        print("::group::Running " + filename)
        result = subprocess.run(
            [sys.executable, "ipynb2py.py", filename],
            capture_output=True, # Optional: captures stdout/stderr
            text=True            # Optional: decodes output as a string
            )

        # Check if it ran successfully
        if result.returncode == 0:
            print("Script executed successfully!")
            print("Output:", result.stdout)
            print("::endgroup::")
        else:

            print("Output:\n", result.stdout)
            print("Error:\n", result.stderr)
            print("::endgroup::")
            print("::error::Running "+filename+" has failed!")
            sys.exit(1)
          
else:
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
        print("A socket error occurred. This is most likely due to a timeout in the NASA Horizons connections. We catch this exception here and ignore it.")
        pass
