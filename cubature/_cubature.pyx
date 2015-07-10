# -*- coding: utf-8 -*-
import numpy as np
cimport numpy as np
import cython
from ._cubature cimport (error_norm, integrand, integrand_v, hcubature, pcubature,
        hcubature_v, pcubature_v)

cdef class Integrand:
    cdef object f, args, kwargs
    cdef unsigned int ndim, fdim

    def __cinit__(self, object f, unsigned ndim, unsigned fdim, *args,
            **kwargs):
        self.f = f
        self.ndim = ndim
        self.fdim = fdim
        self.args = args
        self.kwargs = kwargs

    def __init__(self, *args, **kwargs):
        if not callable(self.f):
            raise ValueError('first argument not callable')
    
    cdef int _call(self, const double *x, double *fval) except -1:
        cdef double [:] _x = <double [:self.ndim]>x 
        cdef double [:] _f = <double [:self.fdim]>fval 
        cdef int error

        try:
            np.asarray(_f)[:] = self.f(np.asarray(_x), *self.args, **self.kwargs)
            error = 0
        except Exception as e:
            error = -1
            raise e
        return error

    cdef int _vcall(self, unsigned npts, const double *x, double *fval) except -1:
        cdef double [:, :] _x = <double [:npts, :self.ndim]>x 
        cdef double [:, :] _f = <double [:npts, :self.fdim]>fval 
        cdef int error

        try:
            np.asarray(_f)[:] = self.f(np.asarray(_x), *self.args, **self.kwargs)
            error = 0
        except Exception as e:
            error = -1
            raise e
        return error

    def call(self, xval):
        cdef double [:] _xval = np.array(xval, dtype=np.float64)
        cdef double [:] fval = np.zeros((self.fdim,), dtype=np.float64)
        err = self._call(&_xval[0], &fval[0])
        if err != 0:
            raise RuntimeError('error while calling call()')
        return np.asarray(fval)

    def vcall(self, xval):
        cdef double [:, :] _xval = np.array(xval, dtype=np.float64)
        cdef unsigned npts = _xval.shape[0]
        cdef double [:, :] fval = np.zeros((npts, self.fdim,), dtype=np.float64)
        err = self._vcall(npts, &_xval[0,0], &fval[0,0])
        if err != 0:
            raise RuntimeError('error while calling vcall()')
        return np.asarray(fval)

cdef int integrand_wrapper(unsigned int ndim, double *x, void *fdata, 
        unsigned int fdim, double *fval) except -1:
    wrapped = <Integrand>fdata;
    return wrapped._call(x, fval) 

cdef int integrand_wrapper_v(unsigned int ndim, unsigned int npts, double *x, 
        void *fdata, unsigned int fdim, double *fval) except -1:
    wrapped = <Integrand>fdata;
    return wrapped._vcall(npts, x, fval) 

@cython.boundscheck(False)
@cython.wraparound(False)
def cubature(callable, unsigned ndim, unsigned fdim, xmin, xmax, str method, 
        double abserr, double relerr, int norm, unsigned maxEval):

    cdef double [:] _xmin = np.array(xmin, dtype=np.float64)
    cdef double [:] _xmax = np.array(xmax, dtype=np.float64)

    cdef double [:] val = np.empty((fdim,), dtype=np.float64)
    cdef double [:] err = np.empty((fdim,), dtype=np.float64)

    wrapper = Integrand(callable, ndim, fdim) 

    if method == 'hcubature_v':
        error = hcubature_v(fdim, <integrand_v>integrand_wrapper_v, 
                <void *> wrapper, ndim, &_xmin[0], &_xmax[0], maxEval, abserr, 
                relerr, <error_norm> norm, &val[0], &err[0])

    elif method == 'hcubature':
        error = hcubature(fdim, <integrand>integrand_wrapper, <void *> wrapper,
                ndim, &_xmin[0], &_xmax[0], maxEval, abserr, relerr, 
                <error_norm> norm, &val[0], &err[0])

    elif method == 'pcubature_v':
        error = pcubature_v(fdim, <integrand_v>integrand_wrapper_v, 
                <void *> wrapper, ndim, &_xmin[0], &_xmax[0], maxEval, abserr,
                relerr, <error_norm> norm, &val[0], &err[0])

    elif method == 'pcubature':
        error = pcubature(fdim, <integrand>integrand_wrapper, <void *> wrapper,
                ndim, &_xmin[0], &_xmax[0], maxEval, abserr, relerr, 
                <error_norm> norm, &val[0], &err[0])

    else:
        raise ValueError('unknown integration method `{:s}`'.format(method))

    if error != 0:
        raise RuntimeError('integration failed')

    return np.asarray(val), np.asarray(err)
