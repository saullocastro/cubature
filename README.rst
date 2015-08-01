.. image:: https://raw.github.com/saullocastro/cubature/master/cubature_logo.png
    :align: center

========
Cubature
========

.. contents::

What is Cubature?
-----------------

It is a numerical integration technique.  From
`MathWorld <http://mathworld.wolfram.com/Cubature.html>`_, 
Ueberhuber (1997, p. 71) and Krommer and Ueberhuber 
(1998, pp. 49 and 155-165) use the word "quadrature" to mean numerical
computation of a univariate integral, and "cubature" to mean numerical
computation of a multiple integral.

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

To install in the ``site-packages`` directory and make it importable from
anywhere:

.. code::
   
    python setup.py install

If you are changing the ``_cubature.pyx`` file, you must have Cython
installed in order to create a new ``_cubature.c`` file. The ``setup.py``
script will automatically try to use the Cython compiler first.

If you want to build only a local ``_cubature.pyd`` file, go to
``./cubature`` and type:

.. code::
   
    python setup.py build_ext -i

Running the tests
-----------------

The Python wrapper has been proven using test integrands from the C
package and some additional testing functions from Genz. The integrands 
were implemented in Cython and verified with Mathematica. 

After building cubature, run the unit tests with the ```pytest``` package in 
the package directory. Be aware that this takes several minutes.

Examples
--------

Some examples are given in `./examples <https://github.com/saullocastro/cubature/tree/master/examples>`_.

Fork me!
--------

You are welcome to fork this repository and modify it in whatever way you
want. It will also be nice if you could send a pull request here in case
you think your modifications is valuable for another person.

License
-------

This wrapper follows the GNU-GPL license terms discribed in the
`C Package <https://github.com/saullocastro/cubature/tree/master/cubature/cpackage/COPYING>`_.
