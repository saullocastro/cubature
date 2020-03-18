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


**Note**: Coverage not taking into account Cython / C files


    
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

Installation from source code
...

You must have Cython installed. Then do::

   python setup.py install 

or (usually in Linux)::

   python3 setup.py install


Installation from pip repository
...

Just do::

   pip install cubature

or (usually in Linux)::

   pip3 install cubature


Running the tests
-----------------

To run the tests you will have to download the source code. After installing as
explained above, go to the source code root folder and run::

    py.test .

The Python wrapper has been proven using test integrands from the C
package and some additional testing functions from Genz. The integrands
were implemented in Cython and verified with Mathematica.


Citing this Python wrapper for Cubature
---------------------------------------

We kindly ask you to cite this Python library properly. Also, it would be
helpful if you could cite the papers where this methods has been applied as
well.

Castro, S.G.P.; Loukianov, A.; et al. "Python wrapper for Cubature: adaptive multidimensional integration". DOI:10.5281/zenodo.2541552. Version 0.14.3, 2020.



Citing Papers using this Python wrapper for Cubature
.....................................................

Used to integrate tangent stiffness matrices in computational solid mechanics

Castro, S.G.P. et al. "Evaluation of non-linear buckling loads of geometrically imperfect composite cylinders and cones with the Ritz method". Composite Structures, Vol. 122, 284-299, 2015.

Castro, S.G.P. et al. "A semi-analytical approach for linear and non-linear analysis of unstiffened laminated composite cylinders and cones under axial, torsion and pressure loads". Thin-Walled Structures, Vol. 90, 61-73, 2015.

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

This wrapper follows the GNU-GPL license terms of Steven G. Johnson described in the `C Package <https://github.com/saullocastro/cubature/tree/master/cubature/cpackage/COPYING>`_.
