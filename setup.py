import os
from subprocess import Popen
from distutils.core import setup
from distutils.extension import Extension

#TODO I am sure there is a more elegant way to do this
cwd = os.getcwd()
os.chdir( os.path.join( os.getcwd(), 'cubature' ) )
p = Popen('python setup.py build_ext -i clean', shell=True)
p.wait()
os.chdir( cwd )
#

setup(
cmdclass = {},
name = 'Cubature',
version = '0.11.0',
description = 'Numerical integration technique',
packages = ['cubature'],
package_data = {'cubature': ['__init__.py', '_cubature.pyd',
                             '_cubature.pyx',
                             'cubature.py', 'test_cubature.py']},
author = 'Saullo G. P. Castro',
author_email = 'saullogiovani@gmail.com',
license = 'GNU-GPL',
keywords =  'numerical integration, integration, '
            'multi-dimensional integration',
url = 'https://github.com/saullocastro/cubature',
long_description = open('README.rst').read(),
)
