import numpy as np
cimport numpy as np
import cython
from cpython cimport tuple, bool, array
from _cubature cimport (error_norm, integrand, integrand_v, hcubature,
                        pcubature, hcubature_v, pcubature_v)

cdef class Integrand:
    cdef object f, data
    cdef unsigned int ndim, fdim

    def __cinit__(self, object f, unsigned ndim, unsigned fdim):
        self.f = f
        self.ndim = ndim
        self.fdim = fdim

    def __init__(self, *args, **kwargs):
        if not callable(self.f):
            raise ValueError('first argument not callable')
    
    cdef int call(self, const double *x, double *fval):
        cdef double [:] _x = <double [:self.ndim]>x 
        cdef double [:] _f = <double [:self.fdim]>fval 
        cdef int error

        try:
            _f = self.f(_x)
            error = 0
        except Exception:
            error = -1

        return error

    cdef int vcall(self, unsigned npts, const double *x, double *fval):
        cdef double [:, :] _x = <double [:npts, :self.ndim]>x 
        cdef double [:, :] _f = <double [:npts, :self.fdim]>fval 
        cdef int error

        try:
            _f = self.f(_x)
            error = 0
        except Exception:
            error = -1

        return error

cdef int integrand_wrapper(unsigned int ndim, double *x, void *fdata, 
        unsigned int fdim, double *fval):
    wrapped = <Integrand>fdata;
    return wrapped.call(x, fval) 

cdef int integrand_wrapper_v(unsigned int ndim, unsigned int npts, double *x, 
        void *fdata, unsigned int fdim, double *fval):
    wrapped = <Integrand>fdata;
    return wrapped.vcall(npts, x, fval) 

@cython.boundscheck(False)
@cython.wraparound(False)
def _cubature(callable,
              unsigned ndim,
              unsigned fdim,
              double [:] xmin,
              double [:] xmax,
              str method,
              double abserr, 
              double relerr, 
              error_norm norm,
              unsigned maxEval,
              double [:] val,
              double [:] err
             ):

    wrapper = Integrand(callable, ndim, fdim) 

    if method == 'hcubature_v':
        ans =  hcubature_v(fdim, <integrand_v>integrand_wrapper_v,
                               <void *> wrapper,
                               ndim,
                               &xmin[0],
                               &xmax[0],
                               maxEval,
                               abserr,
                               relerr,
                               <error_norm> norm,
                               &val[0],
                               &err[0])
    elif method == 'hcubature':
        ans =  hcubature(fdim, <integrand>integrand_wrapper,
                               <void *> wrapper,
                               ndim,
                               &xmin,
                               &xmax,
                               maxEval,
                               abserr,
                               relerr,
                               <error_norm> norm,
                               &val,
                               &err)
    elif method == 'pcubature_v':
        ans =  pcubature_v(fdim, <integrand_v>integrand_wrapper_v,
                               <void *> wrapper,
                               ndim,
                               &xmin[0],
                               &xmax[0],
                               maxEval,
                               abserr,
                               relerr,
                               <error_norm> norm,
                               &val[0],
                               &err[0])
    elif method == 'pcubature':
        ans =  pcubature(fdim, <integrand>integrand_wrapper,
                               <void *> wrapper,
                               ndim,
                               &xmin[0],
                               &xmax[0],
                               maxEval,
                               abserr,
                               relerr,
                               <error_norm> norm,
                               &val[0],
                               &err[0])
    else:
        raise ValueError('unknown integration method')

    return ans
