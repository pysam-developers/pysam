#! /usr/bin/python

'''pysam --- a Python package for reading, manipulating, and writing
genomic data sets.

pysam is a lightweight wrapper of the HTSlib API and provides facilities
to read and write SAM/BAM/CRAM/VCF/BCF/BED/GFF/GTF/FASTA/FASTQ files
as well as access to the command-line functionality of samtools and bcftools.
The module supports compression and random access through indexing.

This module provides a low-level wrapper around HTSlib's C API using Cython
and a high-level API for convenient access to the data within standard genomic
file formats.
'''

import collections
import glob
import logging
import os
import platform
import re
import subprocess
import sys
import sysconfig
import distutils
from contextlib import contextmanager
from setuptools import setup, Command
from setuptools.command.sdist import sdist
from setuptools.extension import Extension

try:
    from setuptools.errors import LinkError
except ImportError:
    from distutils.errors import LinkError

try:
    from Cython.Distutils import build_ext
except ImportError:
    from setuptools.command.build_ext import build_ext

try:
    import cython  # noqa
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

IS_PYTHON3 = sys.version_info.major >= 3
IS_DARWIN = platform.system() == 'Darwin'

log = logging.getLogger('pysam')


@contextmanager
def changedir(path):
    save_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(save_dir)


def run_configure(option):
    sys.stdout.flush()
    is_msys2 = os.environ.get('MSYSTEM') in ('MINGW32', 'MINGW64', 'UCRT64', 'CLANG64', 'CLANG32', 'CLANGARM64')
    
    if is_msys2:
        command = "sh ./configure"
    else:
        command = "./configure"
    
    try:
        retcode = subprocess.call(
            " ".join((command, option)),
            shell=True)
        if retcode != 0:
            return False
        else:
            return True
    except OSError as e:
        return False


def run_make(targets):
    sys.stdout.flush()
    subprocess.check_call([os.environ.get("MAKE", "make")] + targets)


def run_make_print_config():
    stdout = subprocess.check_output(["make", "-s", "print-config"])
    if IS_PYTHON3:
        stdout = stdout.decode("ascii")

    make_print_config = {}
    for line in stdout.splitlines():
        if "=" in line:
            row = line.split("=")
            if len(row) == 2:
                make_print_config.update(
                    {row[0].strip(): row[1].strip()})
    return make_print_config

def run_nm_defined_symbols(objfile):

    is_windows = platform.system() == 'Windows' and 'MSYS' not in platform.release()
    
    command = ["nm", "-g", "-P", objfile]
    
    stdout = subprocess.check_output(command)
    stdout = stdout.decode("utf-8")

    symbols = set()
    for line in stdout.splitlines():
        (sym, symtype) = line.split()[:2]
        if symtype not in "UFNWw":
            if not is_windows:  # This check is adjusted to correctly handle MSYS2 environment.
                # On macOS (Darwin), all symbols have a leading underscore
                if platform.system() == 'Darwin':
                    symbols.add(sym.lstrip('_'))
                else:
                    # Ignore symbols such as _edata (present in all shared objects)
                    if sym[0] not in "_$.@": symbols.add(sym)

    return symbols

# This function emulates the way distutils combines settings from sysconfig,
# environment variables, and the extension being built. It returns a dictionary
# representing the usual set of variables, suitable for writing to a generated
# file or for running configure (provided the returned LIBS is ignored).
def build_config_dict(ext):
    def env(var):
        return [os.environ[var]] if var in os.environ else []

    def sc(var):
        value = sysconfig.get_config_var(var)
        return [value] if value is not None else []

    def optionise(option, valuelist):
        def quote(s): return "'"+s+"'" if " " in s else s
        return list(quote(option+v) for v in valuelist)

    def kvtuples(pairlist):
        def appendoptvalue(t): return t[0] if t[1] is None else t[0]+"="+t[1]
        return map(appendoptvalue, pairlist)

    # For CC, select the first of these that is set
    cc = (env('CC') + sc('CC') + ['gcc'])[0]

    # distutils ignores sysconfig for CPPFLAGS
    cppflags = " ".join(env('CPPFLAGS') + optionise('-I', ext.include_dirs) +
                        optionise('-D', kvtuples(ext.define_macros)) +
                        optionise('-U', ext.undef_macros))

    cflags = " ".join(sc('CFLAGS') + env('CFLAGS') + sc('CCSHARED') +
                      ext.extra_compile_args)

    # distutils actually includes $CPPFLAGS here too, but that's weird and
    # unnecessary for us as we know the output LDFLAGS will be used correctly
    ldflags = " ".join(sc('LDFLAGS') + env('LDFLAGS') + env('CFLAGS') +
                       optionise('-L', ext.library_dirs) +
                       ext.extra_link_args)

    # ext.libraries is computed (incorporating $LIBS etc) during configure
    libs = " ".join(optionise('-l', ext.libraries))

    return { 'CC': cc, 'CPPFLAGS': cppflags, 'CFLAGS': cflags,
             'LDFLAGS': ldflags, 'LIBS': libs }


def write_configvars_header(filename, ext, prefix):
    config = build_config_dict(ext)
    if prefix != 'HTS':
        config['HTSDIR'] = '(unused)'
        config['CURSES_LIB'] = '(unused)'

    log.info("creating %s for '%s' extension", filename, ext.name)
    with open(filename, "w") as outf:
        for var, value in config.items():
            outf.write('#define {}_{} "{}"\n'.format(prefix, var, value))

@contextmanager
def set_compiler_envvars():
    tmp_vars = []
    for var in ['CC', 'CFLAGS', 'LDFLAGS']:
        if var in os.environ:
            if var == 'CFLAGS' and 'CCSHARED' in sysconfig.get_config_vars():
                os.environ[var] += ' ' + sysconfig.get_config_var('CCSHARED')
            print("# pysam: (env) {}={}".format(var, os.environ[var]))
        elif var in sysconfig.get_config_vars():
            value = sysconfig.get_config_var(var)
            if var == 'CFLAGS' and 'CCSHARED' in sysconfig.get_config_vars():
                value += ' ' + sysconfig.get_config_var('CCSHARED')
            print("# pysam: (sysconfig) {}={}".format(var, value))
            os.environ[var] = value
            tmp_vars += [var]

    try:
        yield
    finally:
        for var in tmp_vars:
            del os.environ[var]


def configure_library(library_dir, env_options=None, options=[]):

    configure_script = os.path.join(library_dir, "configure")

    on_rtd = os.environ.get("READTHEDOCS") == "True"
    # RTD has no bzip2 development libraries installed:
    if on_rtd:
        env_options = "--disable-bz2"

    if not os.path.exists(configure_script):
        raise ValueError(
            "configure script {} does not exist".format(configure_script))

    with changedir(library_dir), set_compiler_envvars():
        if env_options is not None:
            if run_configure(env_options):
                return env_options

        for option in options:
            if run_configure(option):
                return option

    return None


def get_pysam_version():
    sys.path.insert(0, "pysam")
    import version
    return version.__version__

if platform.system() == 'Windows':
    def find_msys2_root():

            try:
                df_output = subprocess.check_output(['df', '-a'], text=True)
            except Exception as e:
                print(f"Failed to execute 'df -a': {e}")
                return None

            for line in df_output.splitlines():
                if line.endswith(' /'):
                    # Extract the filesystem, which should be the MSYS2 root
                    filesystem = line.split()[0]
                    return filesystem
            return None

# Override sdist command to ensure Cythonized *.c files are included.
class cythonize_sdist(sdist):
    # Remove when setuptools (as installed on GH runners) has these options
    if not any(opt[0] == 'owner=' for opt in sdist.user_options):
        sdist.user_options.append(('owner=', 'u', 'Specify owner inside tar'))
    if not any(opt[0] == 'group=' for opt in sdist.user_options):
        sdist.user_options.append(('group=', 'g', 'Specify group inside tar'))

    def run(self):
        from Cython.Build import cythonize
        cythonize(self.distribution.ext_modules)
        sdist.run(self)


# Override Cythonised build_ext command to customise macOS shared libraries.

class CyExtension(Extension):
    def __init__(self, *args, **kwargs):
        self._init_func = kwargs.pop("init_func", None)
        self._prebuild_func = kwargs.pop("prebuild_func", None)
        Extension.__init__(self, *args, **kwargs)

    def extend_includes(self, includes):
        self.include_dirs.extend(includes)

    def extend_macros(self, macros):
        self.define_macros.extend(macros)

    def extend_extra_objects(self, objs):
        self.extra_objects.extend(objs)


class cy_build_ext(build_ext):
    def check_ext_symbol_conflicts(self):
        """Checks for symbols defined in multiple extension modules,
        which can lead to crashes due to incorrect functions being invoked.
        Avoid by adding an appropriate #define to import/pysam.h or in
        unusual cases adding another rewrite rule to devtools/import.py.
        """
        symbols = dict()
        for ext in self.distribution.ext_modules:
            for sym in run_nm_defined_symbols(self.get_ext_fullpath(ext.name)):
                symbols.setdefault(sym, []).append(ext.name.lstrip('pysam.'))

        errors = 0
        for (sym, objs) in symbols.items():
            if (len(objs) > 1):
                log.error("conflicting symbol (%s): %s", " ".join(objs), sym)
                #errors += 1

        if errors > 0: raise LinkError("symbols defined in multiple extensions")

    def run(self):
        if sys.platform == 'darwin':
            ldshared = os.environ.get('LDSHARED', sysconfig.get_config_var('LDSHARED'))
            os.environ['LDSHARED'] = ldshared.replace('-bundle', '')

        build_ext.run(self)
        try:
            if HTSLIB_MODE != 'separate':
                self.check_ext_symbol_conflicts()
        except OSError as e:
            log.warning("skipping symbol collision check (invoking nm failed: %s)", e)
        except subprocess.CalledProcessError:
            log.warning("skipping symbol collision check (invoking nm failed)")

    def build_extension(self, ext):

        if isinstance(ext, CyExtension) and ext._init_func:
            ext._init_func(ext)

        if not self.inplace:
            ext.library_dirs.append(os.path.join(self.build_lib, "pysam"))

        if sys.platform == 'darwin':
            # The idea is to give shared libraries an install name of the form
            # `@rpath/<library-name.so>`, and to set the rpath equal to
            # @loader_path. This will allow Python packages to find the library
            # in the expected place, while still giving enough flexibility to
            # external applications to link against the library.
            relative_module_path = ext.name.replace(".", os.sep) + (sysconfig.get_config_var('EXT_SUFFIX') or sysconfig.get_config_var('SO'))
            library_path = os.path.join(
                "@rpath", os.path.basename(relative_module_path)
            )

            if not ext.extra_link_args:
                ext.extra_link_args = []
            ext.extra_link_args += ['-dynamiclib',
                                    '-rpath', '@loader_path',
                                    '-Wl,-headerpad_max_install_names',
                                    '-Wl,-install_name,%s' % library_path,
                                    '-Wl,-x']
        else:
            if not ext.extra_link_args:
            	if platform.system() == "Windows":
                    ext.extra_link_args = ['-shared']
            	else:
                    ext.extra_link_args = []

            ext.extra_link_args += ['-Wl,-rpath,$ORIGIN']

        if isinstance(ext, CyExtension) and ext._prebuild_func:
            ext._prebuild_func(ext, self.force)

        build_ext.build_extension(self, ext)


class clean_ext(Command):
    description = "clean up Cython temporary files"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        objs = glob.glob(os.path.join("pysam", "libc*.c"))
        if objs:
            log.info("removing 'pysam/libc*.c' (%s Cython objects)", len(objs))
        for obj in objs:
            os.remove(obj)

        headers = (glob.glob(os.path.join("htslib",   "*config*.h")) +
                   glob.glob(os.path.join("samtools", "*config*.h")) +
                   glob.glob(os.path.join("bcftools", "*config*.h")))
        if headers:
            log.info("removing '*/*config*.h' (%s generated headers)", len(headers))
        for header in headers:
            os.remove(header)

        objects = (glob.glob(os.path.join("htslib", "*.[oa]")) +
                   glob.glob(os.path.join("htslib", "cram", "*.o")) +
                   glob.glob(os.path.join("htslib", "htscodecs", "htscodecs", "*.o")))
        if objects:
            log.info("removing 'htslib/**/*.o' and libhts.a (%s objects)", len(objects))
        for obj in objects:
            os.remove(obj)


# How to link against HTSLIB
# shared:   build shared chtslib from builtin htslib code.
# external: use shared libhts.so compiled outside of
#           pysam
# separate: use included htslib and include in each extension
#           module. No dependencies between modules and works with
#           setup.py install, but wasteful in terms of memory and
#           compilation time. Fallback if shared module compilation
#           fails.

HTSLIB_MODE = os.environ.get("HTSLIB_MODE", "shared")
HTSLIB_LIBRARY_DIR = os.environ.get("HTSLIB_LIBRARY_DIR", None)
HTSLIB_INCLUDE_DIR = os.environ.get("HTSLIB_INCLUDE_DIR", None)
HTSLIB_CONFIGURE_OPTIONS = os.environ.get("HTSLIB_CONFIGURE_OPTIONS", None)
HTSLIB_SOURCE = None

package_list = ['pysam',
                'pysam.include',
                'pysam.include.samtools',
                'pysam.include.bcftools']
package_dirs = {'pysam': 'pysam',
                'pysam.include.samtools': 'samtools',
                'pysam.include.bcftools': 'bcftools'}

# list of config files that will be automatically generated should
# they not already exist or be created by configure scripts in the
# subpackages.
config_headers = ["samtools/config.h",
                  "bcftools/config.h"]

# If cython is available, the pysam will be built using cython from
# the .pyx files. If no cython is available, the C-files included in the
# distribution will be used.
if HAVE_CYTHON:
    print("# pysam: cython is available - using cythonize if necessary")
    source_pattern = "pysam/libc%s.pyx"
else:
    print("# pysam: no cython available - using pre-compiled C")
    source_pattern = "pysam/libc%s.c"

# Exit if there are no pre-compiled files and no cython available
fn = source_pattern % "htslib"
if not os.path.exists(fn):
    raise ValueError(
        "no cython installed, but can not find {}."
        "Make sure that cython is installed when building "
        "from the repository"
        .format(fn))

print("# pysam: htslib mode is {}".format(HTSLIB_MODE))
print("# pysam: HTSLIB_CONFIGURE_OPTIONS={}".format(
    HTSLIB_CONFIGURE_OPTIONS))
htslib_configure_options = None

if HTSLIB_MODE in ['shared', 'separate']:
    package_list += ['pysam.include.htslib',
                     'pysam.include.htslib.htslib']
    package_dirs.update({'pysam.include.htslib':'htslib'})

    htslib_configure_options = configure_library(
        "htslib",
        HTSLIB_CONFIGURE_OPTIONS,
        ["--enable-libcurl",
         "--disable-libcurl"])

    HTSLIB_SOURCE = "builtin"
    print("# pysam: htslib configure options: {}".format(
        str(htslib_configure_options)))

    config_headers += ["htslib/config.h"]
    if htslib_configure_options is None:
        # create empty config.h file
        with open("htslib/config.h", "w") as outf:
            outf.write(
                "/* empty config.h created by pysam */\n")
            outf.write(
                "/* conservative compilation options */\n")

    with changedir("htslib"):
        htslib_make_options = run_make_print_config()

    for key, value in htslib_make_options.items():
        print("# pysam: htslib_config {}={}".format(key, value))

if platform.system() == 'Windows':
    external_htslib_libraries = ['deflate.dll', 'z', 'bz2.dll', 'crypto.dll', 'lzma.dll', 'regex.dll', 'tre.dll', 'intl.dll', 'iconv.dll', 'ws2_32', 'curl.dll', 'crypt32']
else:
    external_htslib_libraries = ['z']
    if "LIBS" in htslib_make_options:
        external_htslib_libraries.extend(
            [re.sub("^-l", "", x) for x in htslib_make_options["LIBS"].split(" ") if x.strip()])

if HTSLIB_LIBRARY_DIR:
    # linking against a shared, externally installed htslib version,
    # no sources or built libhts.a required for htslib
    htslib_objects = []
    separate_htslib_objects = []
    chtslib_sources = []
    htslib_library_dirs = [HTSLIB_LIBRARY_DIR]
    htslib_include_dirs = [HTSLIB_INCLUDE_DIR]
    if platform.system() == 'Windows':
        external_htslib_libraries = ['hts', 'regex.dll', 'z.dll', 'curl.dll', 'deflate.dll', 'crypto', 'lzma', 'bz2']
    else:
        external_htslib_libraries = ['z', 'hts']
elif HTSLIB_MODE == 'separate':
    # add to each pysam component a separately compiled
    # htslib
    htslib_objects = ['htslib/libhts.a']
    separate_htslib_objects = ['htslib/libhts.a']
    htslib_library_dirs = []
    htslib_include_dirs = ['htslib']
elif HTSLIB_MODE == 'shared':
    # link each pysam component against the same
    # htslib built from sources included in the pysam
    # package.

    # Link with the object files rather than the final htslib/libhts.a, to ensure that
    # all object files are pulled into the link, even those not used by htslib itself.
    htslib_objects = [os.path.join("htslib", x)
                      for x in htslib_make_options["LIBHTS_OBJS"].split(" ")]
    separate_htslib_objects = []

    htslib_library_dirs = ["."] # when using setup.py develop?
    htslib_include_dirs = ['htslib']
else:
    raise ValueError("unknown HTSLIB value '%s'" % HTSLIB_MODE)

# build config.py
with open(os.path.join("pysam", "config.py"), "w") as outf:
    outf.write('HTSLIB = "{}"\n'.format(HTSLIB_SOURCE))
    config_values = collections.defaultdict(int)

    if HTSLIB_SOURCE == "builtin":
        with open(os.path.join("htslib", "config.h")) as inf:
            for line in inf:
                if line.startswith("#define"):
                    key, value = re.match(
                        r"#define (\S+)\s+(\S+)", line).groups()
                    config_values[key] = value
            for key in ["ENABLE_GCS",
                        "ENABLE_PLUGINS",
                        "ENABLE_S3",
                        "HAVE_COMMONCRYPTO",
                        "HAVE_HMAC",
                        "HAVE_LIBBZ2",
                        "HAVE_LIBCURL",
                        "HAVE_LIBDEFLATE",
                        "HAVE_LIBLZMA",
                        "HAVE_MMAP"]:
                outf.write("{} = {}\n".format(key, config_values[key]))
                print("# pysam: config_option: {}={}".format(key, config_values[key]))

# create empty config.h files if they have not been created automatically
# or created by the user:
for fn in config_headers:
    if not os.path.exists(fn):
        with open(fn, "w") as outf:
            outf.write(
                "/* empty config.h created by pysam */\n")
            outf.write(
                "/* conservative compilation options */\n")

#######################################################
include_os = []
os_c_files = []
# for python 3.4, see for example
# http://stackoverflow.com/questions/25587039/
# error-compiling-rpy2-on-python3-4-due-to-werror-
# declaration-after-statement
extra_compile_args = [
    "-Wno-unused",
    "-Wno-strict-prototypes",
    "-Wno-sign-compare",
    "-Wno-error=declaration-after-statement"]

define_macros = []

if platform.system() == "Windows":
    suffix = distutils.sysconfig.get_config_var('EXT_SUFFIX') or distutils.sysconfig.get_config_var('SO')
else:
    suffix = sysconfig.get_config_var('EXT_SUFFIX') or sysconfig.get_config_var('SO')

if platform.system() == 'Windows':
    internal_htslib_libraries = [
        ":libchtslib{}".format(suffix)]
    internal_samtools_libraries = [
        ":libcsamtools{}".format(suffix),
        ":libcbcftools{}".format(suffix),
        ]
    internal_pysamutil_libraries = [
        ":libcutils{}".format(suffix)]
else:
    internal_htslib_libraries = [
        os.path.splitext("chtslib{}".format(suffix))[0]]
    internal_samtools_libraries = [
        os.path.splitext("csamtools{}".format(suffix))[0],
        os.path.splitext("cbcftools{}".format(suffix))[0],
        ]
    internal_pysamutil_libraries = [
        os.path.splitext("cutils{}".format(suffix))[0]]
    external_htslib_objects = []

if platform.system() == 'Windows':
    msys2_root = find_msys2_root()
    if msys2_root:
        external_htslib_objects = []
        external_htslib_objects.append('htslib/libhts.a')
        msystem = subprocess.check_output(['echo', '%MSYSTEM%'], shell=True, text=True).strip()

        if msystem in ['MINGW32', 'MINGW64', 'UCRT64', 'CLANG64', 'CLANG32', 'CLANGARM64']:
            base_path = f"{msys2_root}/{msystem.lower()}/lib"
        else:
            # Default to MINGW64 or handle error/exception as needed
            base_path = f"{msys2_root}/mingw64/lib"

        for lib in external_htslib_libraries:
            if lib != "ws2_32":
                external_htslib_objects.append(f"{base_path}/lib{lib}.a")
        
        # Append ws2_32 library path separately as an exception
        external_htslib_objects.append(f"{base_path}/libws2_32.a")

    else:
        print("MSYS2 root not found.")

if platform.system() == "Windows":
    libraries_for_pysam_module = internal_htslib_libraries + internal_pysamutil_libraries
else:
    libraries_for_pysam_module = external_htslib_libraries + internal_htslib_libraries + internal_pysamutil_libraries

# Order of modules matters in order to make sure that dependencies are resolved.
# The structures of dependencies is as follows:
# libchtslib: htslib utility functions and htslib itself if builtin is set.
# libcsamtools: samtools code (builtin)
# libcbcftools: bcftools code (builtin)
# libcutils: General utility functions, depends on all of the above
# libcXXX (pysam module): depends on libchtslib and libcutils

# The list below uses the union of include_dirs and library_dirs for
# reasons of simplicity.

def prebuild_libchtslib(ext, force):
    if HTSLIB_MODE not in ['shared', 'separate']: return

    write_configvars_header("htslib/config_vars.h", ext, "HTS")

    if force or not os.path.exists("htslib/libhts.a"):
        log.info("building 'libhts.a'")
        with changedir("htslib"):
            # TODO Eventually by running configure here, we can set these
            # extra flags for configure instead of hacking on ALL_CPPFLAGS.
            args = " ".join(ext.extra_compile_args)
            run_make(["ALL_CPPFLAGS=-I. " + args + " $(CPPFLAGS)", "lib-static"])
    else:
        log.warning("skipping 'libhts.a' (already built)")


def prebuild_libcsamtools(ext, force):
    write_configvars_header("samtools/samtools_config_vars.h", ext, "SAMTOOLS")

if platform.system() == 'Windows':
    modules = [
        dict(name="pysam.libchtslib",
             prebuild_func=prebuild_libchtslib,
             sources=[source_pattern % "htslib", "pysam/htslib_util.c"] + os_c_files,
             extra_objects=htslib_objects,
             extra_link_args=["-Wl,--export-all-symbols"]),
        dict(name="pysam.libcsamtools",
             prebuild_func=prebuild_libcsamtools,
             sources=[source_pattern % "samtools"] + glob.glob(os.path.join("samtools", "*.pysam.c")) +
             [os.path.join("samtools", "lz4", "lz4.c")] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=internal_htslib_libraries,
             extra_link_args=["-Wl,--export-all-symbols"]),
        dict(name="pysam.libcbcftools",
             sources=[source_pattern % "bcftools"] + glob.glob(os.path.join("bcftools", "*.pysam.c")) + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=internal_htslib_libraries,
             extra_link_args=["-Wl,--export-all-symbols"]),
        dict(name="pysam.libcutils",
             sources=[source_pattern % "utils", "pysam/pysam_util.c"] + os_c_files,
             extra_objects=separate_htslib_objects,
             #libraries=external_htslib_libraries + internal_htslib_libraries + internal_samtools_libraries),
             libraries=internal_htslib_libraries + internal_samtools_libraries),
        dict(name="pysam.libcalignmentfile",
             sources=[source_pattern % "alignmentfile"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcsamfile",
             sources=[source_pattern % "samfile"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcalignedsegment",
             sources=[source_pattern % "alignedsegment"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libctabix",
             sources=[source_pattern % "tabix"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcfaidx",
             sources=[source_pattern % "faidx"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcbcf",
             sources=[source_pattern % "bcf"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcbgzf",
             sources=[source_pattern % "bgzf"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libctabixproxies",
             sources=[source_pattern % "tabixproxies"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcvcf",
             sources=[source_pattern % "vcf"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
    ]
else:
    modules = [
        dict(name="pysam.libchtslib",
             prebuild_func=prebuild_libchtslib,
             sources=[source_pattern % "htslib", "pysam/htslib_util.c"] + os_c_files,
             extra_objects=htslib_objects,
             libraries=external_htslib_libraries),
        dict(name="pysam.libcsamtools",
             prebuild_func=prebuild_libcsamtools,
             sources=[source_pattern % "samtools"] + glob.glob(os.path.join("samtools", "*.pysam.c")) +
             [os.path.join("samtools", "lz4", "lz4.c")] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=external_htslib_libraries + internal_htslib_libraries),
        dict(name="pysam.libcbcftools",
             sources=[source_pattern % "bcftools"] + glob.glob(os.path.join("bcftools", "*.pysam.c")) + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=external_htslib_libraries + internal_htslib_libraries),
        dict(name="pysam.libcutils",
             sources=[source_pattern % "utils", "pysam/pysam_util.c"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=external_htslib_libraries + internal_htslib_libraries + internal_samtools_libraries),
        dict(name="pysam.libcalignmentfile",
             sources=[source_pattern % "alignmentfile"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcsamfile",
             sources=[source_pattern % "samfile"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcalignedsegment",
             sources=[source_pattern % "alignedsegment"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libctabix",
             sources=[source_pattern % "tabix"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcfaidx",
             sources=[source_pattern % "faidx"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcbcf",
             sources=[source_pattern % "bcf"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcbgzf",
             sources=[source_pattern % "bgzf"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libctabixproxies",
             sources=[source_pattern % "tabixproxies"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
        dict(name="pysam.libcvcf",
             sources=[source_pattern % "vcf"] + os_c_files,
             extra_objects=separate_htslib_objects,
             libraries=libraries_for_pysam_module),
    ]

common_options = dict(
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros,
    # for out-of-tree compilation, use absolute paths
    library_dirs=[os.path.abspath(x) for x in ["pysam"] + htslib_library_dirs],
    include_dirs=[os.path.abspath(x) for x in ["pysam"] + htslib_include_dirs + \
                  ["samtools", "samtools/lz4", "bcftools", "."] + include_os])

if platform.system() == "Windows":
    common_options["extra_objects"] = external_htslib_objects

# add common options (in python >3.5, could use n = {**a, **b}
for module in modules:
    module.update(**common_options)

classifiers = """
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

metadata = {
    'name': "pysam",
    'version': get_pysam_version(),
    'description': "Package for reading, manipulating, and writing genomic data",
    'long_description': __doc__,
    'long_description_content_type': "text/x-rst",
    'author': "Andreas Heger",
    'author_email': "andreas.heger@gmail.com",
    'license': "MIT",
    'platforms': ["POSIX", "UNIX", "MacOS"],
    'classifiers': [_f for _f in classifiers.split("\n") if _f],
    'url': "https://github.com/pysam-developers/pysam",
    'packages': package_list,
    'ext_modules': [CyExtension(**opts) for opts in modules],
    'cmdclass': {'build_ext': cy_build_ext, 'clean_ext': clean_ext, 'sdist': cythonize_sdist},
    'package_dir': package_dirs,
    'package_data': {'': ['*.pxd', '*.h', 'py.typed', '*.pyi'], },
    # do not pack in order to permit linking to csamtools.so
    'zip_safe': False,
}

if __name__ == '__main__':
    dist = setup(**metadata)
