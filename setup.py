import re
import os
import ast

from setuptools import setup, find_packages
import numpy
from distutils.extension import Extension
from Cython.Build import cythonize

_version_re = re.compile(r'__version__\s+=\s+(.*)')

CLASSIFIERS = """\

Development Status :: 5 - Production/Stable
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
Topic :: Scientific/Engineering :: Mathematics
License :: OSI Approved :: BSD License
Programming Language :: Python :: 3.6
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Operating System :: Microsoft :: Windows
Operating System :: Unix

"""

install_requires = [
        "numpy",
        "cython",
        ]

with open('cubature/__init__.py', 'rb') as f:
    version = str(ast.literal_eval(_version_re.search(
        f.read().decode('utf-8')).group(1)))

extensions = [
    Extension('cubature._cubature',
        sources = [
            'cubature/cpackage/hcubature.c',
            'cubature/cpackage/pcubature.c',
            'cubature/get_ptr.c',
            'cubature/_cubature.pyx',
            ],
        include_dirs = [numpy.get_include()],
        language='c',
        ),
    Extension('cubature._test_integrands',
        sources = ['cubature/_test_integrands.pyx'],
        include_dirs = [numpy.get_include()],
        language='c',
        ),
]

ext_modules = cythonize(extensions,
        compiler_directives={'linetrace': True},
        language_level='3',
        )

setup(
    name = 'cubature',
    version = version,
    packages = find_packages(),
    ext_modules = ext_modules,
    include_dirs = [numpy.get_include()],
    author = 'Saullo G. P. Castro and Anton Loukianov',
    author_email = 'saullogiovani@gmail.com',
    license = 'GNU-GPL',
    keywords =  'numerical integration, integration, '
                'multi-dimensional integration',
    url = 'https://github.com/saullocastro/cubature',
    description = 'Numerical integration technique',
    long_description = open('README.rst').read(),
    long_description_content_type = 'text/x-rst',
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    install_requires=install_requires,
)
