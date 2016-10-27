import os
import re
import sys

try:
    from Cython.Distutils import build_ext
except ImportError:
    from setuptools.command.build_ext import build_ext

from distutils.extension import Extension
from distutils.sysconfig import get_config_vars, get_python_lib, get_python_version
from pkg_resources import Distribution


if sys.platform == 'darwin':
    config_vars = get_config_vars()
    config_vars['LDSHARED'] = config_vars['LDSHARED'].replace('-bundle', '')
    config_vars['SHLIB_EXT'] = '.so'
    config_vars['SO'] = '.so'


def is_pip_install():
    if "_" in os.environ and os.environ["_"].endswith("pip"):
        return True
    if "pip-egg-info" in sys.argv:
        return True
    if re.search("/pip-.*-build/", __file__):
        return True
    return False


class CyExtension(Extension):
    def __init__(self, *args, **kwargs):
        self._init_func = kwargs.pop("init_func", None)
        Extension.__init__(self, *args, **kwargs)

    def extend_includes(self, includes):
        self.include_dirs.extend(includes)

    def extend_macros(self, macros):
        self.define_macros.extend(macros)

    def extend_extra_objects(self, objs):
        self.extra_objects.extend(objs)


class cy_build_ext(build_ext):

    def _get_egg_name(self):
        ei_cmd = self.get_finalized_command("egg_info")
        return Distribution(
            None, None, ei_cmd.egg_name, ei_cmd.egg_version, get_python_version(),
            self.distribution.has_ext_modules() and self.plat_name).egg_name()

    def build_extension(self, ext):

        if isinstance(ext, CyExtension) and ext._init_func:
            ext._init_func(ext)

        if not self.inplace:
            ext.library_dirs.append(os.path.join(self.build_lib, "pysam"))

        if sys.platform == 'darwin':

            relative_module_path = ext.name.replace(".", os.sep) + get_config_vars()["SO"]

            if "develop" in sys.argv or "test" in sys.argv:
                # develop-mode and tests use local directory
                pkg_root = os.path.dirname(__file__)
                linker_path = os.path.join(pkg_root, relative_module_path)
            elif "bdist_wheel" in sys.argv or is_pip_install():
                # making a wheel, or pip is secretly involved
                linker_path = os.path.join("@rpath", relative_module_path)
            else:
                # making an egg: `python setup.py install` default behavior
                egg_name = '%s.egg' % self._get_egg_name()
                linker_path = os.path.join("@rpath", egg_name, relative_module_path)

            if not ext.extra_link_args:
                ext.extra_link_args = []
            ext.extra_link_args += ['-dynamiclib',
                                    '-rpath', get_python_lib(),
                                    '-Wl,-headerpad_max_install_names',
                                    '-Wl,-install_name,%s' % linker_path,
                                    '-Wl,-x']
        else:
            if not ext.extra_link_args:
                ext.extra_link_args = []

            ext.extra_link_args += ['-Wl,-rpath,$ORIGIN']
                                    
        build_ext.build_extension(self, ext)
