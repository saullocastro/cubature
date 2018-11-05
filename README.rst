.. image:: https://raw.github.com/saullocastro/cubature/master/cubature_logo.png
    :align: center

========
Cubature
========

|Build status|

.. |Build status| image:: https://travis-ci.org/saullocastro/cubature.svg?branch=master
    :target: https://travis-ci.org/saullocastro/cubature

|Coverage|

.. |Coverage| image:: https://coveralls.io/repos/github/saullocastro/cubature/badge.svg?branch=master
     :target: https://coveralls.io/github/saullocastro/cubature?branch=master

    
.. contents::

What is Cubature?
-----------------

It is a numerical integration technique.  From
`MathWorld <http://mathworld.wolfram.com/Cubature.html>`_,
Ueberhuber (1997, p. 71) and Krommer and Ueberhuber
(1998, pp. 49 and 155-165) use the word "quadrature" to mean numerical
computation of a univariate integral, and "cubature" to mean numerical
computation of a multiple integral.

Documentation
-------------

Please, `see the module documentation here
<http://saullocastro.github.io/cubature/>`_.

Python wrapper for the Cubature package
---------------------------------------

From the `Nanostructures and Computation Wiki at MIT
<http://ab-initio.mit.edu/wiki/index.php/Cubature>`_, `Steven W. Johnson
<http://math.mit.edu/~stevenj/>`_ has written a simple C package for
adaptive multidimensional integration (cubature) of vector-valued
functions over hypercubes and this is a
Python wrapper for the referred C package.

Installation
------------

If you are changing ``_cubature.pyx``, you must have Cython installed in order
to create a new ``_cubature.c`` file (the same is valid for
``_test_integrands.pyx``). The ``setup.py`` script will automatically try to
use the Cython compiler first.

Install to Python's site-package
................................

To install in the ``site-packages`` directory and make it importable from
anywhere:

.. code::

    python setup.py install

Install to a customized site-package
....................................

Windows:

.. code::

    set prefix=anydirectory
    set PYTHONPATH=%anydirectory%\Lib\site-packages;%PYTHONPATH%
    mkdir %anydirectory%
    python setup.py install --prefix="%anydirectory%"

Linux:

.. code::

    export prefix=anydirectory
    export PYTHONPATH=$anydirectory/Lib/site-packages;$PYTHONPATH
    mkdir $anydirectory
    python setup.py install --prefix=$anydirectory


It will create an ``.egg`` file that will go to
``$anydirectory\Lib\site-packages``.  This file can be unzipped to obtain the
importable module, OR, the full path to this ``.egg`` can be added to
``$PYTHONPATH`` (which can be also done inside a script through
sys.path.append().


Build locally
.............

If you want to build it locally (without installing in Python's
``site-packages``) just type:

.. code::

    python setup.py build_ext --inplace clean

Running the tests
-----------------

The Python wrapper has been proven using test integrands from the C
package and some additional testing functions from Genz. The integrands
were implemented in Cython and verified with Mathematica.

After building cubature, run the unit tests with the ```pytest``` package in
the package directory. Be aware that this takes several minutes:

.. code::

    python -m pytest test_cubature.py


Cite Cubature
--------------

We kindly ask you to cite this Python library propertly. Also, it would be
helpful if you could cite the papers where this methods has been applied as
well.

Papers: Used to integrate tangent stiffness matrices in computational solid mechanics
**************************************************************************************

Castro, S.G.P. et al. "Evaluation of non-linear buckling loads of geometrically imperfect
composite cylinders and cones with the Ritz method". Composite Structures, Vol. 122, 284-299, 2015.

Castro, S.G.P. et al. "A semi-analytical approach for linear and non-linear analysis of unstiffened laminated composite cylinders and cones under axial, torsion and pressure loads". Thin-Walled Structures, Vol. 90, 61-73, 2015.

The Python wrapper should be cited as follows
*********************************************

Castro, S.G.P.; Loukianov, A.; et al. "Python wrapper for Cubature: adaptive multidimensional integration". On-line: https://github.com/saullocastro/cubature/releases, Version 0.13.3, 2017.
(check if the version and year are correct)

Examples
--------

Some examples are given in `./examples <https://github.com/saullocastro/cubature/tree/master/examples>`_.

Fork me!
--------

You are welcome to fork this repository and modify it in whatever way you
want. It will also be nice if you could send a pull request here in case
you think your modifications are valuable for another person.

License
-------

This wrapper follows the GNU-GPL license terms discribed in the
`C Package <https://github.com/saullocastro/cubature/tree/master/cubature/cpackage/COPYING>`_.
