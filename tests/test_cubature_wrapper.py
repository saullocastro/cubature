import ctypes
import numpy as np
import pytest

import importlib
cubature_module = importlib.import_module('cubature.cubature')
from cubature import cubature as cub


def test_invalid_xmin_length():
    def f(x):
        return np.array([0.0])
    with pytest.raises(AssertionError):
        cub(f, ndim=2, fdim=1, xmin=[0], xmax=[1, 1])


def test_invalid_vectorized_output_shape():
    def f(x):
        # return wrong shape: not (7, fdim)
        return np.ones((5, 2))
    with pytest.raises(ValueError, match=r"shape=\(:, fdim\)"):
        cub(f, ndim=1, fdim=2, xmin=[0], xmax=[1], vectorized=True)


def test_invalid_scalar_vectorized_output():
    def f(x):
        # return 2D array but fdim=1
        return np.ones((7, 1, 1))
    with pytest.raises(ValueError, match=r"valid array"):
        cub(f, ndim=1, fdim=1, xmin=[0], xmax=[1], vectorized=True)


def test_invalid_adaptive_option():
    def f(x):
        return np.array([0.0])
    with pytest.raises(ValueError, match="unknown combination"):
        cub(f, ndim=1, fdim=1, xmin=[0], xmax=[1], adaptive='z')


def test_raw_callback(monkeypatch):
    called = {}

    def fake_raw(func, ndim, fdim, xmin, xmax, method, abserr, relerr, norm, maxEval, args=(), kwargs=None):
        called['raw'] = True
        return np.array([1.0]), np.array([0.0])

    def fake_standard(*a, **k):
        called['standard'] = True
        return np.array([1.0]), np.array([0.0])

    monkeypatch.setattr(cubature_module, '_cython_cubature_raw_callback', fake_raw)
    monkeypatch.setattr(cubature_module, '_cython_cubature', fake_standard)

    CBTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_void_p)

    def c_func(x):
        return 0.0

    c_cb = CBTYPE(c_func)

    cub(c_cb, ndim=1, fdim=1, xmin=[0], xmax=[1])

    assert 'raw' in called and 'standard' not in called
