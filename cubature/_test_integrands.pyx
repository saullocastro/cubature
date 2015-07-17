# cython: profile=False
import numpy as np
cimport numpy as np
from libc.math cimport cos
from libc.math cimport M_PI as pi

cimport cython
#from cpython.array cimport array
#cdef array double_template = array('d')

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef genz_oscillatory(np.ndarray x,  np.ndarray a, double u):
    return np.cos(2*np.pi*u + np.sum(a*x))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double genz_oscillatory_fast(unsigned int n, double *args):
    # there should be 2*d + 1 parameters
    cdef unsigned int d = (n - 1)//2
    cdef double [:] vargs = <double [:n]>args
    cdef double [:] x = vargs[:d]
    cdef double [:] a = vargs[d:n-1]
    cdef double u = vargs[n-1]

    cdef double val = 0.

    for i in range(d): 
        val += a[i]*x[i]
    
    val = cos(2*pi*u + val)

    return val

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef genz_oscillatory_c(np.ndarray x,  np.ndarray a, double u):
    cdef unsigned int d = x.shape[0]
    cdef unsigned int n = 2*d + 1
    cdef np.ndarray[dtype=np.float64_t, ndim=1] args = np.empty((n,),
            dtype=float)
    args[:d] = x[:d]
    args[d:2*d] = a[:d]
    args[n-1] = u
    return genz_oscillatory_fast(n, <double *>args.data)

cpdef genz_oscillatory_exact(double n, np.ndarray a,
        double u):
    '''calculate the exact integral of the Genz oscillatory function on the
    half interval (0, n) given an array of frequencies a and phase u'''

    cdef unsigned m = a.shape[0]
    return 2**m *np.cos(0.5 * np.sum(a) * n + 2*np.pi*u) * \
            np.prod(np.sin(0.5 * a * n))/np.prod(a)
