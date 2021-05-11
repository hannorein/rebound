# Installation

## Choosing between C and Python

You can use either C or Python when working with REBOUND.
Which programming language you want to use depends on your taste and your specific application. In short: 

- If you want to set up a planetary system, visualize data with matplotlib, and integrate your simulation with one of the built-in integrators then use the Python version. It's quick and easy to use. 
- If you want to run large simulations with millions of particles, develop your own integrator, use OpenGL visualizations, or make use of the distributed tree code, then use the C version. C gives you the best performance and direct access to all the internals.

!!! Note
    All the computationally expensive parts of REBOUND are written in C. So even if you use the Python version, your simulation will run very efficiently.
    If you want to extend REBOUND, for example to include an additional non-gravitational force, you can do that in both C or Python. For complicated force routines, a C implementation of your function would most likely be significantly faster.  


## Installation via pip
If you just want to try out REBOUND or don't plan to modify it in any way, then the easiest way to install the python version of REBOUND is [pip](https://pypi.org) (the Package Installer for Python). Simply type the following command into a terminal:

```bash
pip install rebound
```

If you have trouble installing a package with pip, consider using a [virtual environment](https://docs.python.org/3/tutorial/venv.html).

## Installation via git

We use the [git](https://git-scm.com) as a version control system for REBOUND. 
If you want to use the C version of REBOUND or plan to make any modifications to REBOUND, you can clone the repository to your computer. 
Make sure you have git installed, then type the following command in a terminal:

``` bash
git clone https://github.com/hannorein/rebound
```

This will create a new directory names `rebound/` which contains all the source code, examples, and documentation.
To use the python version of REBOUND, go to the `rebound/` directory, then install REBOUND with 
```bash
pip install -e .
```
You should now be able to import REBOUND from python. 

!!! Info "Installing REBOUND on Windows"
    To best way to use REBOUND on Windows is to install the Windows Subsystem (WSL) for Linux.
    The Ubuntu distribution works well with REBOUND.
    Once you've installed WSL, open a WSL terminal and make sure you have a compiler installed with the command `sudo apt-get install gcc`.
    If you want to use the Python version of REBOUND, also make sure you have a recent version of Python installed: `sudo apt-get install python3 python3-pip`.
    From here on, you can follow the installation instructions above.

## Compiling the C version of REBOUND

### Examples

If you look at any of the examples in the `examples/` subdirectories, you'll see one
`problem.c` file and one `Makefile`. All the REBOUND code itself is in the
`src/` directory. This setup keeps the different projects nicely separated from the shared REBOUND code.
To compile one of the examples, go to the example's directory and type `make`. 
This triggers the following tasks:

1.  The Makefile in the example directory sets up various environment variables. These
    determine settings like compiler optimization flags and which
    libraries are included (see below).
2.  Next, the Makefile in the `src/` directory gets called. This compiles
    the entire REBOUND code into a shared library.
3.  It then creates a symbolic link from the current directory to the
    location of the shared library is located.
4.  Finally, it compiles your own code, the `problem.c` file and links it to the REBOUND shared library.

You can execute your program with `./rebound`. After you edited either the `problem.c` file or any file in the `src/` directory, you can simply type `make` again to recompile your program. 
If you change any of the environment variables, clean the build directory first, by executing `make clean`.

### Your own project

The easiest way to start working on your own problem is to simply copy an example directory that is somewhat similar to what you want to do.
This way, all your project's source and data files will be in one directory, separate from the main REBOUND source files in `src/`. 

Alternatively, you can also install the shared REBOUND library in a global directory (e.g. `/usr/lib/`) and the header file in `/usr/include/`. 
Doing so will allow you (and any other users on your system) to use REBOUND from any directory.
However, doing so requires root access and some knowledge on how Unix systems work.
By simply replicating and modifying one of the examples, you'll avoid these complications.

### Possible issues during compilation

The way we've designed REBOUND should make the compilation process extremely easy.
You do not need to install any additional libraries (although you might want to, see below), and you do not need root access. 
You might nevertheless run into problems. Some of the most common issues are:

-   **Missing compilers.** Make sure you have a C compiler installed. If
    you are using a Mac, install the Xcode package which you can
    download for free on the App Store. Make sure the command line tools 
    are installed.
-   **Missing glfw3 library.** You can compile REBOUND with support for
    real-time OpenGL visualizations. This requires the glfw3 library. If
    you are on a Mac, then the easiest way to install the glfw3 library
    is with homebrew:
    `brew tap homebrew/versions && brew install glfw3`. If you are on
    Linux, you can install it with your package manager, for example
    with `sudo apt-get install libglfw3-dev`. Alternatively, you can
    disable the OpenGL visualization in the Makefile by setting
    `OPENGL=0`. Then, execute `make clean` and try compiling the program
    again. Note that on some systems the `glfw` library is called
    `glfw3` instead. In that case, change `-lglfw` to `-lglfw3` 
    in the file `src/Makefile.defs`.
-   **Issue with march=native.** Some users have reported issues related
    to the compiler flag `-march=native` which tries to optimize the
    code for the native architecture. If you want to have the 
    most optimized code, add the `-march=native` or `-mtune=native` flag
    in the file `src/Makefile.defs`. If you use the python version, you
    can add compiler flags to `setup.py`. This might improve performance
    significantly.
