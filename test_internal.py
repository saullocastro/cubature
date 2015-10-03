from __future__ import division

from cubature._cubature import Integrand
import numpy as np

def test_call_1arg():
    def func(x, a, b, c, kw1=True, kw2=None):
        return 1.

    d = 10
    a = Integrand(func, d, 1, args=(2, 3, 'a'), kwargs={'kw1': 1, 'kw2': 2})
    a.call(np.zeros((2,), dtype=float))

def test_call_narg():
    ndim = 1
    fdim = 10
    def func(x, a, b, c, kw1=True, kw2=None):
        return np.zeros((fdim,))

    a = Integrand(func, ndim, fdim, args=(2, 3, 'a'), kwargs={'kw1': 1, 'kw2': 2})
    a.call(np.zeros((1,), dtype=float))

def test_vcall_1arg():
    def func(x, a, b, c, kw1=True, kw2=None):
        return np.ones((x.shape[0],))

    d = 10
    a = Integrand(func, d, 1, args=(2, 3, 'a'), kwargs={'kw1': 1, 'kw2': 2})
    a.vcall(np.zeros((100, d), dtype=float))

def test_vcall_narg():
    ndim = 1
    fdim = 10
    def func(x, a, b, c, kw1=True, kw2=None):
        return np.zeros((x.shape[0], fdim,))

    a = Integrand(func, ndim, fdim, args=(2, 3, 'a'), kwargs={'kw1': 1, 'kw2': 2})
    a.vcall(np.zeros((1000, 1,), dtype=float))
