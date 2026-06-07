# This file exists solely to configure the dynamic C extension.
# All project metadata lives in pyproject.toml.

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from glob import glob
import platform
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


##### AVX512 where supported
def is_x86_64():
    machine = platform.machine().lower()
    return "x86_64" in machine or "amd64" in machine

class build_ext_avx512(build_ext):
    def build_extensions(self):
        os.makedirs(self.build_temp, exist_ok=True)
        if is_x86_64():
            asm = os.path.join(self.build_temp, "integrator_whfast512.asm_o")
            self.compiler.spawn(["as", "-g", "-o", asm, "src/integrator_whfast512.s"])
            for ext in self.extensions:
                ext.extra_objects = (ext.extra_objects or []) + [asm]

        old_compile = self.compiler._compile

        def custom_compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
            # Clone args to prevent cross-contamination between files
            current_postargs = list(extra_postargs) if extra_postargs else []

            # Check if this specific compilation pass is targeting ARM64
            is_arm64_pass = "-arch arm64" in " ".join(
                cc_args
            ) or "-arch arm64" in " ".join(current_postargs)
            is_x86_pass = "-arch x86_64" in " ".join(
                cc_args
            ) or "-arch x86_64" in " ".join(current_postargs)

            # Handle your assembly file routing
            if src.endswith(".s"):
                if is_arm64_pass and not is_x86_pass:
                    print(f"--- Skipping x86 assembly {src} for ARM64 pass ---")
                    return

                if is_x86_pass:
                    current_postargs.extend(
                            [
                                "-fno-integrated-as",
                                "-Wa,-march=generic64+avx512f+avx512dq+avx512bw+avx512vl",
                                ]
                            )

            return old_compile(
                obj, src, ext, cc_args, current_postargs, pp_opts
            )

        # Swap the compiler's compile hook with our filtered architecture rule
        self.compiler._compile = custom_compile
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
