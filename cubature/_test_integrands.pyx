import numpy as np
cimport numpy as np
from libc.math cimport cos
from libc.math cimport M_PI as pi

cimport cython
from cpython.array cimport array, clone
cdef array double_template = array('d')

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
