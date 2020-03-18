#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False
#cython: infer_types=False

cimport numpy as np
import numpy as np
import cython
from cpython.ref cimport PyObject
from ._cubature cimport (error_norm, integrand, integrand_v, hcubature, pcubature,
        hcubature_v, pcubature_v)

cdef extern from "get_ptr.h":
    void *get_ctypes_function_pointer(PyObject *obj)

cdef class Integrand:
    cdef object f, args, kwargs
    cdef unsigned int ndim, fdim

    def __cinit__(self, object f, unsigned ndim, unsigned fdim, object args,
            object kwargs):
        self.f = f
        self.ndim = ndim
        self.fdim = fdim
        self.args = args
        self.kwargs = kwargs

    def __init__(self, *args, **kwargs):
        if not callable(self.f):
            raise ValueError('first argument not callable')

    def __str__(self):
        s = 'Integrand(f = {!r}, ndim = {!r}, fdim = {!r}, args = {!r}, kwargs = {!r})'\
             .format(self.f, self.ndim, self.fdim, self.args, self.kwargs)
        return s

    cdef int _call(self, const double *x, double *fval) except -1:
        cdef double [:] _x = <double [:self.ndim]>x
        cdef double [:] _f = <double [:self.fdim]>fval
        cdef int error
        cdef unsigned i

        try:
            tmp = self.f(np.asarray(_x), *self.args, **self.kwargs)
            if self.fdim == 1:
                _f[0] = tmp
            else:
                for i in range(self.fdim):
                    _f[i] = tmp[i]
            error = 0
        except Exception as e:
            error = -1
            raise e
        return error

    cdef int _vcall(self, unsigned npts, const double *x, double *fval) except -1:
        cdef double [:, :] _x = <double [:npts, :self.ndim]>x
        cdef double [:, :] _f = <double [:npts, :self.fdim]>fval
        cdef int error
        cdef unsigned i,j

        try:
            tmp = self.f(np.asarray(_x), *self.args, **self.kwargs)
            if self.fdim == 1:
                for i in range(npts):
                    _f[i] = tmp[i]
            else:
                for i in range(npts):
                    for j in range(self.fdim):
                        _f[i,j] = tmp[i,j]
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
        unsigned int fdim, double *fval):
    cdef Integrand wrapped = <Integrand>fdata;
    return wrapped._call(x, fval)


cdef int integrand_wrapper_v(unsigned int ndim, unsigned int npts, double *x,
        void *fdata, unsigned int fdim, double *fval):
    cdef Integrand wrapped = <Integrand>fdata;
    return wrapped._vcall(npts, x, fval)


def cubature(callable, unsigned ndim, unsigned fdim, xmin, xmax, str method,
        double abserr, double relerr, int norm, unsigned maxEval, args=(),
        kwargs={}):

    cdef double [:] _xmin = np.array(xmin, dtype=np.float64)
    cdef double [:] _xmax = np.array(xmax, dtype=np.float64)

    cdef double [:] val = np.empty((fdim,), dtype=np.float64)
    cdef double [:] err = np.empty((fdim,), dtype=np.float64)

    wrapper = Integrand(callable, ndim, fdim, args, kwargs)

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
        raise ValueError('unknown integration method `{!s}`'.format(method))

    if error != 0:
        raise RuntimeError('integration failed')

    return np.asarray(val), np.asarray(err)

def cubature_raw_callback(callable, unsigned ndim, unsigned fdim, xmin, xmax, str method,
        double abserr, double relerr, int norm, unsigned maxEval, args=(),
        kwargs={}):

    cdef double [:] _xmin = np.array(xmin, dtype=np.float64)
    cdef double [:] _xmax = np.array(xmax, dtype=np.float64)

    cdef double [:] val = np.empty((fdim,), dtype=np.float64)
    cdef double [:] err = np.empty((fdim,), dtype=np.float64)

    if method == 'hcubature_v':
        error = hcubature_v(fdim, <integrand_v>get_ctypes_function_pointer(<PyObject *>callable),
                NULL, ndim, &_xmin[0], &_xmax[0], maxEval, abserr,
                relerr, <error_norm> norm, &val[0], &err[0])

    elif method == 'hcubature':
        error = hcubature(fdim, <integrand>get_ctypes_function_pointer(<PyObject *>callable), NULL,
                ndim, &_xmin[0], &_xmax[0], maxEval, abserr, relerr,
                <error_norm> norm, &val[0], &err[0])

    elif method == 'pcubature_v':
        error = pcubature_v(fdim, <integrand_v>get_ctypes_function_pointer(<PyObject *>callable),
                NULL, ndim, &_xmin[0], &_xmax[0], maxEval, abserr,
                relerr, <error_norm> norm, &val[0], &err[0])

    elif method == 'pcubature':
        error = pcubature(fdim, <integrand>get_ctypes_function_pointer(<PyObject *>callable), NULL,
                ndim, &_xmin[0], &_xmax[0], maxEval, abserr, relerr,
                <error_norm> norm, &val[0], &err[0])

    else:
        raise ValueError('unknown integration method `{!s}`'.format(method))

    if error != 0:
        raise RuntimeError('integration failed')

    return np.asarray(val), np.asarray(err)
