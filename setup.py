# This file exists solely to configure the dynamic C extension.
# All project metadata lives in pyproject.toml.

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from glob import glob
import platform
import os
import subprocess
import sys
import sysconfig

##### Git hash
try:
    ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()
    ghash_arg = f"-DGITHASH={ghash}"
except Exception:
    ghash_arg = "-DGITHASH=cf99043eb72d820c5b131b1b543d9254660fd44d" #GITHASHAUTOUPDATE

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


##### AVX512 where supported
def avx512_supported():
    machine = platform.machine().lower()
    is_x86_64 = "x86_64" in machine or "amd64" in machine
    return is_x86_64

class build_ext_avx512(build_ext):
    def build_extensions(self):
        os.makedirs(self.build_temp, exist_ok=True)
        if sys.platform == "darwin":
            if hasattr(self.compiler, 'linker_so'):
                for i, flag in enumerate(self.compiler.linker_so):
                    if flag == '-bundle':
                        self.compiler.linker_so[i] = '-shared'
        if avx512_supported():
            if sys.platform == "win32":
                orig_spawn = self.compiler.spawn
                def new_spawn(cmd, **kwargs):
                    if cmd and "cl.exe" in cmd[0].lower():
                        cmd[0] = "clang-cl"
                    return orig_spawn(cmd, **kwargs)
                self.compiler.spawn = new_spawn
                asm = os.path.join(self.build_temp, "integrator_whfast512_asm.obj")
                cmd = ["clang", "-c", "src/integrator_whfast512.s", "-o", asm, "-target", "x86_64-pc-windows-msvc"]
                try:
                    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                except subprocess.CalledProcessError as e:
                    print(f"Assembly compilation failed: {e.stderr.decode()}", file=sys.stderr)
                    raise
            else:
                asm = os.path.join(self.build_temp, "integrator_whfast512_asm.o")
                self.compiler.spawn(["as", "--noexecstack", "-g", "-o", asm, "src/integrator_whfast512.s"])
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
