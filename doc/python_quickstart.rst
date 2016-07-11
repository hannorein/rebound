Quick User Guide (Python)
=========================

Installation
------------

Installing REBOUND is very easy. Rebound does not depend on any libraries. However, you need to have python (version 2 or 3) and a C compiler installed. Most likely, you already have those on your system.  

If you don't, and aren't sure how to go about getting them, it is probably easiest to install either the Enthought or Anaconda python distributions (which are free and come with typically used libraries and an easy-to-use installer).  For the C compiler on Mac, it's probably easiest to install Xcode from the App store.

Note:  REBOUND does not work on Windows, and we currently do not have plans to support it.

Create a virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before you start, you can create a virtual environment to keep your python installation clean. This is recommended, but not necessary. You need to have [virtualenv](https://virtualenv.pypa.io/en/latest/) installed (if you use the Anaconda python distribution, you'll need to instead create a conda environment - see below).

To create a virtual environment, simply open a terminal window, go to the folder where you want REBOUND to be installed (e.g. `/home/username/rebound/`) and enter::

    virtualenv venv

This creates a virtual environment named `venv`. To activate it, type::

    source venv/bin/activate

If you log out of your terminal or open a new one, you'll need to reactivate the virtual environment with the above command (if the environment is active you'll see its name in parentheses before the command prompt).

If you use the Anaconda distribution, the above likely will not work. To create a conda environment run the following command::

    conda create -n venv pip
    source activate venv


Standard python installation using pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now you can install REBOUND using pip. All you have to do is type the following into a terminal window (if you created a virtual environment above you should first run the appropriate activate command depending on whether you used virtualenv or conda):

    pip install rebound

The setup script will download the latest version of REBOUND from [PyPI](https://pypi.python.org/pypi) (the Python Package Index), compile the C code in the background and place all the files in their correct location. No other libraries are needed to start working with WHFast and REBOUND, but you might want to install numpy and matplotlib to be able to post-process your data and make plots. For analysis tools (and to run FourierSpectrum.ipynb) you might also want scipy.  Installing those libraries is also very easy (but may need a few minutes).  Depending on whether you use virtualenv or conda, use

    pip install numpy matplotlib scipy
    
or
    
    conda install numpy matplotlib scipy
    

That's all there is to do!


Installing the development version directly from github
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of using pip to install the latest version of REBOUND from PyPI, you can also install REBOUND directly from github. This is useful if you want to make changes to REBOUND.

First clone the repository (you might need to create a github account)::

    git clone git@github.com:hannorein/rebound.git

If you already know that you might be contributing something to REBOUND, you can also first fork the repository on github and then clone your own repository.

Next, create a virtual environment with the same commands as above::

    virtualenv venv && source venv/bin/activate

Now you can install rebound from source using::

    pip install -e .

You can modify the python module files in the directory `rebound/` and you'll see the changed the next time you run a python script (no need to reinstall the REBOUND package every time).

If you install REBOUND directly from github, you can also run it without python. Have a look at the README file in the main directory and at the examples in the `examples/` directory for examples in C. These are much more diverse than the python examples (e.g. allow you to use a tree code for gravity calculations, use other boundary conditions, etc).

Python Notebook
^^^^^^^^^^^^^^^


The tutorials in python_tutorials were written using iPython/Jupyter notebooks. You can view them directly on GitHub. If you want to edit them or create your own notebook, you'll need to install iPython (make sure to activate the virtual environment first if you created one)::

    pip install "ipython[notebook]"

or::

    conda install ipython-notebook
    
You can then open iPython notebooks in your browser by typing::

    ipython notebook
    
To create a new notebook select from the dropdown menu on the top right the item that says 'New'. Now you can interactively follow the commands in the tutorials or run your own!


First REBOUND simulation
------------------------

To run your first REBOUND simulation, just start python with the command::

    python

Then, import the rebound module::

    import rebound

create a new simualtion::

    sim = rebound.Simulation()

Now you can add as many particles to REBOUND as you want::

    sim.add(m=1.0)
    sim.add(m=1.0e-3, a=1.0)

Above, we added a star with mass 1 and a planet with mass 0.001 at 1 AU. By default REBOUND uses units in which G=1. Next you can start integrating your particles forward in time::

    sim.integrate(1000.)

Now, the time has advanced to t=1000. You can print out the particle positions with::

    sim.status()

For more information, have a look at the python examples, which act as tutorials. You can also read the documentation for the REBOUND module and the REBOUND C code to get a better understanding of what is going on behind the scenes.

Upgrading REBOUND
-----------------

REBOUND is actively being expanded and improved, so it's worthwhile to periodically update it.

If you installed REBOUND with::

    pip install rebound

then simply::

    pip install rebound --upgrade

If this does not work, you have an old version of pip.  You can either upgrade pip (probably best!), or simply::

    pip uninstall rebound
    pip install rebound

If you cloned the git repository, i.e., have a `rebound` folder on your file system, merge the new changes (see here_ if you're not sure how to do this) and::

    pip install -e .

.. _here: https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging
