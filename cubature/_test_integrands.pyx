# cython: profile=False

import numpy as np
cimport numpy as np
from libc.math cimport cos, sin, exp, tgamma
from libc.math cimport M_PI as pi

cimport cython
from cpython.array cimport array, clone
cdef array double_template = array('d')

cdef double k2sqrtpi = 1.12837916709551257390 

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double genz_oscillatory(double [:] x,  double [:] a, double u):
    cdef double sum = 0.
    cdef unsigned int i
    for i in range(x.shape[0]):
        sum += x[i]*a[i]
    return np.cos(2*np.pi*u + sum)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double genz_oscillatory_fast(unsigned int n, double *args):
    # there should be 2*d + 1 parameters
    cdef unsigned int d = (n - 1)//2
    cdef double val = 0.

    cdef unsigned int i
    for i in range(d): 
        val += args[i]*args[d+i]
    
    val = cos(2*pi*args[n-1] + val)

    return val

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double genz_oscillatory_c(double [:] x,  double [:] a, double u):
    cdef unsigned int d = x.shape[0]
    #cdef unsigned int n = 2*d + 1
    
    #cdef array[double] args = clone(double_template, n, 0)
    #for i in range(d):
    #    args[i] = x[i]
    #    args[d+i] = a[i]
    #args[n-1] = u
    #return genz_oscillatory_fast(n, args.data.as_doubles)

    cdef double val = 0.
    cdef unsigned int i
    for i in range(d): 
        val += x[i]*a[i]
    
    val = cos(2*pi*u + val)
    return val

cpdef genz_oscillatory_exact(double n, np.ndarray a,
        double u):
    '''calculate the exact integral of the Genz oscillatory function on the
    half interval (0, n) given an array of frequencies a and phase u'''

    cdef unsigned m = a.shape[0]
    return 2**m *np.cos(0.5 * np.sum(a) * n + 2*np.pi*u) * \
            np.prod(np.sin(0.5 * a * n))/np.prod(a)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double cubature_zero(double [:] x):
    '''reimplementation of cubature's test function 0'''
    cdef unsigned int d = x.shape[0]

    cdef double val = 1.
    cdef unsigned int i

    for i in range(d):
        val *= cos(x[i])

    return val

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double cubature_zero_exact(double [:] x):
    '''calculate exact integral of cubature's test function zero on the half
    interval (0, xmax[i]) in d dimensions'''
    cdef unsigned int d = x.shape[0]
    cdef unsigned int i
    cdef double val = 1.

    for i in range(d):
        val *= sin(x[i])

    return val

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double cubature_one(double [:] x):
    '''calculate cubature test function one, a rescaled gaussian'''
    cdef unsigned int d = x.shape[0]
    cdef unsigned int i
    cdef double sum = 0., z, scale = 1.
    
    for i in range(d):
        if x[i] > 0.:
            z = (1 - x[i])/x[i] # rescale to unit interval
            sum += z*z
            scale *= k2sqrtpi / (x[i] * x[i])
        else:
            scale = 0.

    return scale * exp(-sum)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double cubature_one_exact(double [:] x):
    return 1.

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double cubature_two(double [:] x, double radius):
    cdef unsigned int d = x.shape[0]
    cdef double val = 0.

    for i in range(d):
        val += x[i]*x[i]

    return 1. if val < radius*radius else 0.

cdef double nsphere_surface_area(unsigned int d, double radius):
    d = d+1
    return d*pow(pi, 0.5*d)/tgamma(0.5*d + 1) * pow(radius, d)

cpdef double cubature_two_exact(unsigned int d, double radius):
    return nsphere_surface_area(d, radius)/(d + 1)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double cubature_three(double [:] x):
    cdef unsigned int d = x.shape[0]
    cdef double prod = 1.
    cdef unsigned int i

    for i in range(d):
        prod *= 2*x[i]
    return prod

