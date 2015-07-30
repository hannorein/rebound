REBOUND - An open-source multi-purpose N-body code
==================================================

.. image:: http://img.shields.io/badge/rebound-v2.2.1-green.svg?style=flat
.. image:: http://img.shields.io/badge/license-GPL-green.svg?style=flat :target: https://github.com/hannorein/rebound/blob/master/LICENSE
.. image:: http://img.shields.io/travis/hannorein/rebound/master.svg?style=flat :target: https://travis-ci.org/hannorein/rebound/
.. image:: http://img.shields.io/badge/arXiv-1110.4876-orange.svg?style=flat :target: http://arxiv.org/abs/1110.4876
.. image:: http://img.shields.io/badge/arXiv-1409.4779-orange.svg?style=flat :target: http://arxiv.org/abs/1409.4779
.. image:: http://img.shields.io/badge/arXiv-1506.01084-orange.svg?style=flat :target: http://arxiv.org/abs/1506.01084


NEW VERSION
-----------
Welcome to REBOUND version 2! We made many changes to the code. Most importanly, REBOUND is now thread-safe and does not use global variables anymore. All the variables that were previously global, are now contained in the `reb_simulation` structure. This has many advantages, for example, you can run separate simulations in parallel from within one process. We also made it possible to choose all modules at runtime (compared to the selection in the `Makefile` that was used before). This is much more in line with standard UNIX coding practice and does not severely impact performance (it might even help making REBOUND a tiny bit faster). This makes REBOUND a fully functional shared library. We added a prefix to all public functions and struct definitions: `reb_`.

There are still some features that haven't been fully ported. Most importantly, the MPI parallelization and the SWEEP collision detection routine. 

The best way to get and idea of the changes we made is to look at some of the example problems. If you have trouble using the new version or find a bug, please submit an issue or a pull request on github. 

We're working on a unified way of documenting the code. Stay tuned (or even better, help us!).

-------------------



