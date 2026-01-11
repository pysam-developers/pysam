import os
import re
import sys
import sysconfig

try:
    from Cython.Distutils import build_ext
except ImportError:
    from setuptools.command.build_ext import build_ext

from setuptools.extension import Extension


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
        super().__init__(*args, **kwargs)

    def extend_includes(self, includes):
        self.include_dirs.extend(includes)

    def extend_macros(self, macros):
        self.define_macros.extend(macros)

    def extend_extra_objects(self, objs):
        self.extra_objects.extend(objs)


class cy_build_ext(build_ext):

    def run(self):
        if sys.platform == 'darwin':
            ldshared = os.environ.get('LDSHARED', sysconfig.get_config_var('LDSHARED'))
            os.environ['LDSHARED'] = ldshared.replace('-bundle', '')

        super().run()

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
            relative_module_path = ext.name.replace(".", os.sep) + sysconfig.get_config_var("EXT_SUFFIX")
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
                ext.extra_link_args = []

            ext.extra_link_args += ['-Wl,-rpath,$ORIGIN']
                                    
        super().build_extension(ext)
