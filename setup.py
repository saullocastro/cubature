import numpy
from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
    ext_modules = [Extension('cubature._cubature', ['cubature/_cubature.pyx'])]
    cmdclass = {'build_ext': build_ext}
except ImportError:
    ext_modules = [Extension('cubature._cubature', ['cubature/_cubature.c'])]
    cmdclass = {}
    pass

setup(
cmdclass = cmdclass,
ext_modules = ext_modules,
include_dirs = [numpy.get_include()],
name = 'Cubature',
version = '0.11.0',
description = 'Numerical integration technique',
packages = ['cubature'],
author = 'Saullo G. P. Castro',
author_email = 'saullogiovani@gmail.com',
license = 'GNU-GPL',
keywords =  'numerical integration, integration, '
            'multi-dimensional integration',
url = 'https://github.com/saullocastro/cubature',
long_description = open('README.rst').read(),
)
