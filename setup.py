import re
import ast
from setuptools import setup, Extension
import numpy

_version_re = re.compile(r'__version__\s+=\s+(.*)')

CLASSIFIERS = """\

Development Status :: 5 - Production/Stable
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
Topic :: Scientific/Engineering :: Mathematics
License :: OSI Approved :: BSD License
Operating System :: Microsoft :: Windows
Programming Language :: Python :: 3.5
Programming Language :: Python :: 3.6
Programming Language :: Python :: 3.7
Operating System :: Unix

"""

install_requires = [
        "numpy",
        "cython",
        ]

with open('cubature/__init__.py', 'rb') as f:
    version = str(ast.literal_eval(_version_re.search(
        f.read().decode('utf-8')).group(1)))

try:
    from Cython.Distutils import build_ext
    ext_modules = [
        Extension('cubature._cubature',
            sources = [
                'cubature/_cubature.pyx',
                'cubature/cpackage/hcubature.c',
                'cubature/cpackage/pcubature.c',
                'cubature/get_ptr.c'],
            include_dirs = [numpy.get_include()]),
        Extension('cubature._test_integrands',
            sources = ['cubature/_test_integrands.pyx'],
            include_dirs = [numpy.get_include()]),
    ]
    cmdclass = {'build_ext': build_ext}
except ImportError:
    ext_modules = [
            Extension('cubature._cubature', ['cubature/_cubature.c']),
            Extension('cubature._test_integrands', ['cubature/_test_integrands.c'])
            ]
    cmdclass = {}

setup(
    name = 'cubature',
    version = version,
    cmdclass = cmdclass,
    packages = ['cubature'],
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
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    install_requires=install_requires,
)
