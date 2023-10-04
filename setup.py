try:
    from setuptools import setup, Extension
except ImportError:
    # Legacy distutils import. No longer available on Pyton > 3.12
    from distutils.core import setup, Extension
from codecs import open
import os
import sys

import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

# Try to get git hash
try:
    import subprocess
    ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii")
    ghash_arg = "-DGITHASH="+ghash.strip()
except:
    ghash_arg = "-DGITHASH=e9d5c3e37cbf3f235fc6f6bd6bd763076d4773b0" #GITHASHAUTOUPDATE

extra_link_args=[]
if sys.platform == 'darwin':
    config_vars = sysconfig.get_config_vars()
    config_vars['LDSHARED'] = config_vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args=['-Wl,-install_name,@rpath/librebound'+suffix]
if sys.platform == 'win32':
    extra_compile_args=[ghash_arg, '-DLIBREBOUND', '-D_GNU_SOURCE']
else:
    # Default compile args
    extra_compile_args=['-fstrict-aliasing', '-O3','-std=c99','-Wno-unknown-pragmas', ghash_arg, '-DLIBREBOUND', '-D_GNU_SOURCE', '-fPIC']

# Option to disable FMA in CLANG. 
FFP_CONTRACT_OFF = os.environ.get("FFP_CONTRACT_OFF", None)
if FFP_CONTRACT_OFF:
    extra_compile_args.append('-ffp-contract=off')

# Option to enable AVX512 enabled integrators (WHFast512). 
AVX512 = os.environ.get("AVX512", None)
if AVX512:
    extra_compile_args.append('-march=native')
    extra_compile_args.append('-DAVX512')
    
libreboundmodule = Extension('librebound',
                    sources = [ 'src/rebound.c',
                                'src/integrator_ias15.c',
                                'src/integrator_whfast.c',
                                'src/integrator_whfast512.c',
                                'src/integrator_saba.c',
                                'src/integrator_mercurius.c',
                                'src/integrator_eos.c',
                                'src/integrator_leapfrog.c',
                                'src/integrator_bs.c',
                                'src/integrator_janus.c',
                                'src/integrator_sei.c',
                                'src/integrator.c',
                                'src/gravity.c',
                                'src/boundary.c',
                                'src/display.c',
                                'src/collision.c',
                                'src/tools.c',
                                'src/fmemopen.c',
                                'src/rotations.c',
                                'src/derivatives.c',
                                'src/tree.c',
                                'src/particle.c',
                                'src/binarydiff.c',
                                'src/output.c',
                                'src/input.c',
                                'src/simulationarchive.c',
                                'src/transformations.c',
                                ],
                    include_dirs = ['src'],
                    define_macros=[ ('LIBREBOUND', None) ],
                    extra_link_args=extra_link_args,
                    extra_compile_args=extra_compile_args,
                    )

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='rebound',
    version='3.28.3',
    description='An open-source multi-purpose N-body code',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/hannorein/rebound/',
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
