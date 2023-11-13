# Installation

![type:video](https://www.youtube.com/embed/_7Y3YLKyxWA)

## Choosing between C and Python

You can use either C or Python when working with REBOUND.
Which programming language you want to use depends on your preference and your specific application. In short: 

- If you want to set up a planetary system, visualize data with matplotlib, and integrate your simulation with one of the built-in integrators then use the Python version. It's quick and easy to use. 
- If you want to run large simulations with millions of particles, develop your own integrator, use the distributed tree code with MPI, OpenMP parallelization, or OpenGL visualization, then use the C version. C gives you the best performance and direct access to all the REBOUND internals.

!!! Note
    All the computationally expensive parts of REBOUND are written in C. So even if you use the Python version, your simulation will run very efficiently.
    If you want to extend REBOUND, for example to include an additional non-gravitational force, you can do that in both C or Python. However, for complicated force routines, a C implementation of your function would most likely be significantly faster.  


## Installation via pip
!!! info inline end "Python Wheels"
    Starting with REBOUND version 3.28, we provide Python Wheels for REBOUND. 
    This makes installing REBOUND easier on a wide variety of systems. 
    For optimal performance, you can compile REBOUND yourself with optimizations flags that specifically target your system.


If you just want to try out REBOUND or don't plan to modify it in any way, then the easiest way to install the python version of REBOUND is [pip](https://pypi.org) (the Package Installer for Python). Simply type the following command into a terminal:

```bash
pip install rebound
```

If you have trouble installing a package with pip, consider using a [virtual environment](https://docs.python.org/3/tutorial/venv.html).
Also, make sure your version of pip is not too old. You can update pip with pip itself:
```bash
pip install --upgrade pip
```

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

## Compiling the C version of REBOUND

### Examples

If you look at any of the examples in the `examples/` sub-directories, you'll see one
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
    location of the shared library in the `src/` directory. (On Windows
    the  Makefile simply copies the shared library instead of making a symbolic link)
4.  Finally, it compiles your own code, the `problem.c` file and links it to the REBOUND shared library.

You can execute your program with `./rebound` (or `rebound.exe` on Windows). 
After you edit either the `problem.c` file or any file in the `src/` directory, you can simply type `make` again to recompile your program. 
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
    are installed. if you are on Windows, make sure you install the 
    compilers that come with Visual studio (`cl.exe`) and have them available
    in your current command prompt (use the *Developer Command Prompt for VS*).
-   **Missing glfw3 library.** You can compile REBOUND with support for
    real-time OpenGL visualizations. This is an optional feature that 
    requires the glfw3 library. If you are on a Mac, then the easiest 
    way to install the glfw3 library is with homebrew:
    `brew tap homebrew/versions && brew install glfw3`. If you are on
    Linux, you can install it with your package manager, for example
    with `sudo apt-get install libglfw3-dev`. Alternatively, you can
    disable the OpenGL visualization in the Makefile by setting
    `OPENGL=0`. Then, execute `make clean` and try compiling the program
    again. Note that on some systems the `glfw` library is called
    `glfw3` instead. In that case, change `-lglfw` to `-lglfw3` 
    in the file `src/Makefile.defs`.
-   **Compiler optimizations.** By default, REBOUND does not use the
    compiler flag `-march=native` which tries to optimize the
    code for the native architecture. If you want to have the 
    most optimized code, add the `-march=native` or `-mtune=native` flag
    in the file `src/Makefile.defs`. If you use the python version, you
    can add compiler flags to `setup.py`. This might improve performance
    significantly.
-   **Floating point contractions.** Some compilers (e.g. clang) optimize code by
    contracting certain floating point operations (e.g. a multiplication
    and an addition become one fused multiply-add instruction). This improves performance but might prevent you from 
    reproducing results exactly. You can turn off fused multiply-add instruction with the
    `-ffp-contract=off` compiler flag. If you use the python version, you can set the
    `FFP_CONTRACT_OFF` environment variable before installing REBOUND.


## Running REBOUND on Windows

There are several ways to run REBOUND on Windows.

### Python
You can install the python version of REBOUND using pip:
```bash
pip install -e .
```
This will download the latest python wheel and install it on your system. 

### Windows Subsystem for Linux (WSL)
You can run the C-version of REBOUND using the Windows Subsystem for Linux (WSL).
You will need `make` and a compiler, such as `gcc`. These can be installed within WSL with the following command:
```bash
sudo apt-install make gcc
```
Then, you can follow the above instructions for Linux. Start by download REBOUND, for example using git:
```bash
git clone https://github.com/hannorein/rebound
```
Then, compile and run a simple C-example with the following commands:
```bash
cd rebound/examples/simplest
make
./rebound
```

### Native Windows Builds
!!! note inline end Note
    The native Windows support for REBOUND is relatively new. Several features are currently not supported on native Windows builds: OpenMP, MPI, OpenGL, and AVX512. Please [file a bug report on github](https://github.com/hannorein/rebound/issues) if you require any of these featured or if you encounter any other problems. 

Since version 3.28, you can also run REBOUND natively on Windows. You need to install make and enable the Microsoft Visual Studio compiler. Once you have downloaded the source code of REBOUND, open the Developer Command Prompt for VS or the Windows PowerShell on your system and go to the REBOUND source code. Then, compile and run a simple C-example with the following commands:
```bash
cd examples
cd simplest
make
rebound.exe
```

