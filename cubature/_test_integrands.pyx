# cython: profile=False
import numpy as np
cimport numpy as np

cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef genz_oscillatory(np.ndarray x,  np.ndarray a, float u):
    return np.cos(2*np.pi*u + np.sum(a*x))

cpdef genz_oscillatory_exact(float n, np.ndarray a,
        float u):
    '''calculate the exact integral of the Genz oscillatory function on the
    half interval (0, n) given an array of frequencies a and phase u'''

    cdef unsigned m = a.shape[0]
    return 2**m *np.cos(0.5 * np.sum(a) * n + 2*np.pi*u) * \
            np.prod(np.sin(0.5 * a * n))/np.prod(a)
