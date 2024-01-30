# Visualization

Starting with version 4, REBOUND includes new real-time interactive 3D visualizations. 
The code for these visualization is built-upon the previous OpenGL visualizations that came with the C version of REBOUND. 
However, the new visualization feature has several advantages

- There are zero dependencies. No need to install GLFW or other libraries.
- The visualizations work on Linux, MacOS, Windows, including mobile devices.
- The visualizations work both for the C and python version of REBOUND.
- You can use visualizations for simulation running on remote servers.

This page describes how to use these visualizations and the technology that makes this possible behind the scenes.

## Basic idea

Getting 3D visualization work out of the box on a wide variety of platforms is difficult. 
To circumvent most of the issues, REBOUND uses an un-conventional approach: it allows you to use your web browser.
Web browsers have the advantage of being readily available on all operating systems.
And because they implement open standards (HTML, JS, WebAssembly, WebGL) they are unlikely to break compatibility with REBOUND any time soon.

To get data from a simulation to a web browser, **REBOUND comes with its own built-in web server!**
In fact, every simulation can have its own web server. 
After starting a web server, it sits idle in a separate thread until you decide to use it.
When a simulation is deallocated, the server is stopped. 

The following code shows how to start the server:

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_start_server(r, 1234);

    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.start_server(port=1234)
    ```

By default, the server opens port 1234 on your computer. 
Now all you have to do to see the visualization is to open your browser and go to [http://localhost:1234/](http://localhost:1234/) or [http://127.0.0.1:1234/](http://127.0.0.1:1234/).

When you open the page the REBOUND web server accepts your request and serves you a `rebound.html` file which includes all the code required to visualize a simulation using WebGL. 
The cool thing is, the visualization code is just REBOUND itself, compiled to WebAssembly using emscripten. 
So the visualization that you see in your web browser is exactly the same as the one you see when compiling REBOUND with the `OPENGL=1` option but without all the hassles associated with using OpenGL/glut/GLFW libraries.

You can generate (compile) a `rebound.html` file yourself. The code for that is in the directory `web_client/`. 
But to do that you would need to download and install emscripten. 
To help you out when connection to the REBOUND web server, REBOUND first looks for a `rebound.html` file in the current directly. 
If it doesn't find it, then it downloads a pre-compiled file from GitHub.
This should work seamlessly in the background, so you might not even notice.

Now, after REBOUND served the visualization code to your web browser, it needs to be able to send simulation data to the browser.
This is done by packing up the simulation as binary data in the form of a Simulationarchive. 
This data is then sent to your browser via HTTP whenever your browser requests a new frame for the visualization.
This can be up to 60 times per second.
The REBOUND version running in your browser then reconstructs the full simulation using the Simulationarchive data.

You can also send simple commands back from the browser to the server.
For example, by pressing the space bar you can pause and un-pause the integration.
You can press `q` to terminate the integration. 
Commands that only affect the visualization (for example you can press `w` to show/hide orbits) are not sent back to the main simulation on the server.
For all keyboard commands available, press `h`. A help window will show up on screen.


## Widget in Jupyter notebooks
Instead of opening a new browser window for the visualization, you can also include a visualization widget in your Jupyter notebook.

```python
sim = rebound.Simulation()
sim.widget(size=(400,400))
```

This automatically starts the web server and the shows an iframe in your current notebook. 
The iframe is simply showing the contents of http://localhost:1234.
You can connect multiple browser windows to the same simulation. 
So in addition to having the widget directly in your notebook, you can also go to [http://localhost:1234](http://localhost:1234) to see the visualization.



## Multiple simulations
You can visualize multiple simulations at the same time. 
For that to work, each simulation needs to have its own port. 
Make sure you close the server if you want to re-use the port the server is using for another simulation.

=== "C"
    ```c
    reb_simulation_stop_server(r);  // This stops the server.
    reb_simulation_free(r);         // This also stops the server if it's still running.
    ```

=== "Python"
    ```python
    sim.stop_server(port=1234)      # This stops the server.
    del sim                         # This also stops the server if it's still running.
    ```

## Connecting to remote servers

Sometimes you might run a simulation on a remote server or computing cluster. 
When connecting to the remote server using ssh, you can forward the port used for visualization. 
REBOUND uses port 1234 by default, so you might want to enable port forwarding using

```bash
ssh username@remotecomputer -L 1234:localhost:1234
```

You can then connect to the visualization as usual by pointing your browser to http://localhost:1234.

## Disabling web server

Although REBOUND is compiled with the web server capability by default, no web server is started until you call `reb_simulation_start_server()` or `sim.start_sever()`.
You can also disable the web server capability completely if you want.
To do that set `export SERVER=0` in the Makefile.

## Security and resource considerations
The built-in REBOUND web server provides a quick and easy way to visualize simulations. 
It is not intended to be exposed to the public internet because someone might be able to gain access to your computer.

Note that the visualization uses a considerable amount of CPU resources. You might want to disable it (stop the server) if you no longer use it.

!!! info inline end "Future optimizations"
    There are in principle better ways to stream data from a server to a client, for example using WebSockets.
    A future version of REBOUND might optimize the CPU and bandwidth usage.

Depending on your simulation (e.g. for simulations with a large number of particles), the communication between the REBOUND web server and your browser might use a lot of bandwidth and a lot of HTTP requests (up to 60 requests per second). 


