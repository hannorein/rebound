from setuptools import setup, Extension
from codecs import open
from os import path

libreboundmodule = Extension('librebound',
                    sources = [ 'src/librebound.c',
                                'src/integrator_ias15.c',
                                'src/integrator_mikkola.c',
                                'src/integrator_wh.c',
                                'src/integrator_leapfrog.c',
                                'src/integrator_sei.c',
                                'src/integrator_hybrid.c',
                                'src/integrator.c',
                                'src/gravity_direct.c',
                                'src/tools.c',
                                'src/particle.c',
                                'src/output.c',
                                ],
                    include_dirs = ['src'],
                    define_macros=[ ('LIBREBOUND', None), ('COLLISIONS_NONE',None) ],
                    extra_compile_args=['-O3','-std=c99'],
                                    )

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='rebound',
    version='0.1.23',
    description='The N-body integrator REBOUND',
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
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],
    keywords='astronomy astrophysics nbody integrator symplectic wisdom-holman',
    packages=['rebound'],
    install_requires=[],
    ext_modules = [libreboundmodule],
    zip_safe=False)
