try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension
from codecs import open
import os
import sys

import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

extra_link_args=[]
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args=['-Wl,-install_name,@rpath/librebound'+suffix]
    
libreboundmodule = Extension('librebound',
                    sources = [ 'src/rebound.c',
                                'src/integrator_ias15.c',
                                'src/integrator_whfast.c',
                                'src/integrator_wh.c',
                                'src/integrator_leapfrog.c',
                                'src/integrator_sei.c',
                                'src/integrator_hybrid.c',
                                'src/integrator.c',
                                'src/gravity.c',
                                'src/boundary.c',
                                'src/collision.c',
                                'src/tools.c',
                                'src/derivatives.c',
                                'src/tree.c',
                                'src/particle.c',
                                'src/output.c',
                                'src/input.c',
                                ],
                    include_dirs = ['src'],
                    define_macros=[ ('LIBREBOUND', None) ],
                    extra_compile_args=['-fstrict-aliasing', '-O3','-std=c99','-march=native','-Wno-unknown-pragmas', '-DLIBREBOUND', '-D_GNU_SOURCE', '-fPIC'],
                    extra_link_args=extra_link_args,
                    )

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='rebound',
    version='2.16.1',
    description='An open-source multi-purpose N-body code',
    long_description=long_description,
    url='http://github.com/hannorein/rebound',
    author='Hanno Rein',
    author_email='hanno@hanno-rein.de',
    license='GPL',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'Topic :: Scientific/Engineering :: Astronomy',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],
    keywords='astronomy astrophysics nbody integrator symplectic wisdom-holman',
    packages=['rebound'],
    package_data = {'rebound':['rebound.h']},
    install_requires=[],
    tests_require=["numpy","matplotlib"],
    test_suite="rebound.tests",
    ext_modules = [libreboundmodule],
    zip_safe=False)
