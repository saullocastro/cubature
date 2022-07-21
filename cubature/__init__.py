"""
Cubature (:mod:`cubature`)
==========================

.. currentmodule:: cubature

There is one function that embodies all functionalities offered by Cubature.

One of the nice things about Cubature is that you can perform
multi-dimensional integrations at once, which can be achieved by defining:

- `ndim`
- `xmin`
- `xmax`

There are two types of refinement in the integration methods, both necessary
to achieve a desired integration error threshold:

- `h` additional integration points are added
- `p` the order of the integration polynomials is increased

It allows the evaluation of vectorized functions, making it convenient to take
advantage of NumPy's speed (see examples with `vectorized=True`).

See the detailed description below on how to use a vector valued fuction
(function that returs an array) or a scalar function.

.. autofunction:: cubature

More Examples
=============

.. literalinclude:: ../../examples/ex_volumes.py
.. literalinclude:: ../../examples/ex_areas.py



"""
__version__ = '0.14.15'

from .cubature import *

