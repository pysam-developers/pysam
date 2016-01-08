import os
import sys

from distutils.sysconfig import get_config_vars, get_python_lib, get_python_version
from distutils.extension import Extension
from pkg_resources import Distribution
from Cython.Distutils import build_ext

if sys.platform == 'darwin':
    config_vars = get_config_vars()
    #config_vars['SO'] = '.dylib'
    config_vars['LDSHARED'] = config_vars['LDSHARED'].replace('-bundle', '-Wl,-x')


class CyExtension(Extension):
    def __init__(self, name, sources, include_dirs=None, define_macros=None, undef_macros=None, library_dirs=None,
                 libraries=None, runtime_library_dirs=None, extra_objects=None, extra_compile_args=None, extra_link_args=None,
                 export_symbols=None, swig_opts=None, depends=None, language=None, init_func=None, **kw):

        self._init_func = init_func

        Extension.__init__(self, name, sources, include_dirs, define_macros, undef_macros, library_dirs, libraries,
                           runtime_library_dirs, extra_objects, extra_compile_args, extra_link_args, export_symbols,
                           swig_opts, depends, language, **kw)

    def extend_includes(self, includes):
        self.include_dirs.extend(includes)

    def extend_macros(self, macros):
        self.define_macros.extend(macros)

    def extend_extra_objects(self, objs):
        self.extra_objects.extend(objs)


class cy_build_ext(build_ext):
    def build_extension(self, ext):
        if isinstance(ext, CyExtension) and ext._init_func:
            ext._init_func(ext)

        if sys.platform == 'darwin':
            ext.extra_link_args = ext.extra_link_args or []

            egg_name = '%s.egg' % self._get_egg_name()
            name = ext.name
            install_name = '%s/%s/%s%s' % (get_python_lib(), egg_name, name.replace('.', os.sep), config_vars['SO'])

            extra = ['-dynamiclib', '-undefined', 'dynamic_lookup', '-shared',
                     '-Wl,-headerpad_max_install_names', '-Wl,-install_name,%s' % install_name]

            ext.extra_link_args += extra

        build_ext.build_extension(self, ext)

    def _get_egg_name(self):
        ei_cmd = self.get_finalized_command("egg_info")
        return Distribution(None, None, ei_cmd.egg_name, ei_cmd.egg_version, get_python_version(),
                            self.distribution.has_ext_modules() and self.plat_name).egg_name()
