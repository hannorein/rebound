# This file exists solely to configure the dynamic C extension.
# All project metadata lives in pyproject.toml.

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
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
    ghash_arg = "-DGITHASH=53773afe758ed4540dd515968ca9726f851ea5e7" #GITHASHAUTOUPDATE

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
    extra_compile_args = [ "-fstrict-aliasing", "-std=c99", "-Wno-unreachable-code", "-Wno-unknown-pragmas", ghash_arg, "-D_GNU_SOURCE", "-DSERVER", "-fPIC"]
    if os.environ.get("COVERAGE"):
        extra_compile_args += ["-O1", "-fprofile-arcs", "-ftest-coverage", "--coverage"]
        extra_link_args   += ["-fprofile-arcs", "-ftest-coverage", "--coverage"]
    else:
        extra_compile_args.append("-O3")

##### Turn off floating point contractions for bitwise reproducibility
if os.environ.get("FFP_CONTRACT_OFF"):
    extra_compile_args.append("-ffp-contract=off")

AVX512_FILES = {"integrator_asm512.c", "integrator_whfast512.c"}
AVX512_FLAGS = ["-DAVX512", "-mavx512f", "-mavx512dq", "-mavx512cd", "-mavx512bw", "-mavx512vl"]

class build_ext_avx512(build_ext):
    def build_extensions(self):
        base = self.compiler._compile
        def _compile(obj, src, ext, cc_args, postargs, pp_opts):
            if os.path.basename(src) in AVX512_FILES:
                postargs = postargs + AVX512_FLAGS
            return base(obj, src, ext, cc_args, postargs, pp_opts)
        self.compiler._compile = _compile
        os.makedirs(self.build_temp, exist_ok=True)
        asm = os.path.join(self.build_temp, "integrator_asm512.o")
        self.compiler.spawn(["as", "-g", "-o", asm, "src/integrator_asm512.s"])
        for ext in self.extensions:
            ext.extra_objects = (ext.extra_objects or []) + [asm]
        super().build_extensions()

##### C Extension
libreboundmodule = Extension(
    "librebound",
    sources=sorted(glob("src/*.c")),
    include_dirs=["src"],
    extra_link_args=extra_link_args,
    extra_compile_args=extra_compile_args,
)

setup(ext_modules=[libreboundmodule], cmdclass={"build_ext": build_ext_avx512})
