# This file exists solely to configure the dynamic C extension.
# All project metadata lives in pyproject.toml.

from setuptools import setup, Extension
from glob import glob
import os
import sys
import sysconfig

##### Git hash
try:
    import subprocess
    ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()
    ghash_arg = f"-DGITHASH={ghash}"
except Exception:
    ghash_arg = "-DGITHASH=16889b7f9d2381f2aa3ac1571e4a548343ef15ce"  # GITHASHAUTOUPDATE

##### Link args
extra_link_args = []
if sys.platform == "darwin":
    cfg = sysconfig.get_config_vars()
    cfg["LDSHARED"] = cfg["LDSHARED"].replace("-bundle", "-shared")
    suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
    extra_link_args = [f"-Wl,-install_name,@rpath/librebound{suffix}"]

##### Compile args
if sys.platform == "win32":
    extra_compile_args = [ghash_arg, "-DBUILDINGLIBREBOUND", "-D_GNU_SOURCE", "-DSERVER"]
else:
    extra_compile_args = [ "-fstrict-aliasing", "-std=c99", "-Wno-unknown-pragmas", ghash_arg, "-D_GNU_SOURCE", "-DSERVER", "-fPIC"]
    if os.environ.get("COVERAGE"):
        extra_compile_args += ["-O1", "-fprofile-arcs", "-ftest-coverage", "-coverage"]
        extra_link_args   += ["-fprofile-arcs", "-ftest-coverage", "-coverage"]
    else:
        extra_compile_args.append("-O3")

##### Turn off floating point contractions for bitwise reproducibility
if os.environ.get("FFP_CONTRACT_OFF"):
    extra_compile_args.append("-ffp-contract=off")

if os.environ.get("AVX512"):
    extra_compile_args += ["-march=native", "-DAVX512"]

##### C Extension
libreboundmodule = Extension(
    "librebound",
    sources=sorted(glob("src/*.c")),
    include_dirs=["src"],
    extra_link_args=extra_link_args,
    extra_compile_args=extra_compile_args,
)

setup(ext_modules=[libreboundmodule])
