import numpy as np
from ._cubature import cubature as _cython_cubature

__all__ = ['ERROR_INDIVIDUAL', 'ERROR_PAIRED', 'ERROR_L2', 'ERROR_L1', 
        'ERROR_LINF', 'cubature']

ERROR_INDIVIDUAL = 0
ERROR_PAIRED = 1
ERROR_L2 = 2
ERROR_L1 = 3
ERROR_LINF = 4

# maps (adaptive, vectorized) to appropriate function call
# in 
_call_map = {
    ('h', True): 'hcubature_v',
    ('h', False): 'hcubature',
    ('p', True): 'pcubature_v',
    ('p', False): 'pcubature',
    }

def cubature(func, ndim, xmin, xmax, args=tuple(), kwargs=dict(),
             abserr=1.e-8, relerr=1.e-8, norm=ERROR_INDIVIDUAL, maxEval=0,
             adaptive='h', vectorized=False):
    r"""Numerical-integration using cubature technique.
    Parameters
    ----------
    func : callable
        If ``vectorized=False`` the callable must have the form:
            ``f(x_array, *args, **kwargs)``
            - `x_array` in an array containing the `ndim` variables
               being integrated
            - `args` is a tuple containing any other arguments
              required by the function
            - `kwargs` is a dict containing any keyword arguments
              required by the function
            - the function must return a 1-D `np.ndarray` object
            - example::
                  def func(x_array, *args, **kwargs):
                      x, y = x_array
                      return np.array([x**2-y**2, x*y, x*y**2])
        If ``vectorized=True`` the function must have the form:
            ``f(x_array, npt, *args, **kwargs)``
            - `x_array` has ``shape[0]=ndim*npt``
            - `npt` tells the number of points passed to the function
            - `args` is a tuple containing any other arguments
              required by the function
            - `kwargs` is a dict containing any keyword arguments
              required by the function
            - example::
                  def func(x_array, npt, *args, **kwargs):
                      ndim = 2
                      fdim = 3
                      out = np.zeros(fdim*npt)
                      for i in range(npt):
                          x = x_array[i*ndim]
                          y = x_array[i*ndim+1]
                          out[i*fdim] = x**2-y**2
                          out[i*fdim+1] = x*y
                          out[i*fdim+2] = x*y**2
                      return out
        The results from both vectorized and non-vectorized examples
        above should be the same, but the vectorized implementation
        allows more control for the user to parallelize the computation
        of the integration points.
    ndim : integer
        Number dimensions or number of variables being integrated.
    xmin : numpy.ndarray
        A 1-D array carring the minimum integration limit for each
        variable being integrated. It must be have:
        ``xmin.shape[0]=ndim``.
    xmax : numpy.ndarray
        A 1-D array carring the maximum integration limit for each
        variable being integrated. It must be have:
        ``xmax.shape[0]=ndim``.
    args : tuple or list
        Contains the extra arguments required by `func`.
    kwargs : dict-like 
        Contains the extra keyword arguments passed to `func`.
    adaptive : string
        The adaptive scheme used along the adaptive integration.
        'h' means 'h-adaptive', where the domain is partitioned
        'p' means 'p-adaptive', where the order of the integration rule
                                is increased
        The 'p-adaptive' scheme is often better for smoth functions in
        low dimensions.
    epsabs : double
        The maximum number of function evaluations.
    epsrel : double
        The maxmum number of function evaluations.
    norm : integer
        Specifies the norm that is used to measure the error and
        determine convergence properties (irrelevant for
        single-valued functions).
        The `norm` argument takes one of the values:
        - ERROR_INDIVIDUAL: convergence is achieved only when each
            integrand individually satisfies the requested error
            tolerances;
        - ERROR_PAIRED: like ERROR_INDIVIDUAL, except that the
            integrands are grouped into consecutive pairs, with the error
            tolerance applied in a L2 sense to each pair. This option is
            mainly useful for integrating vectors of complex numbers,
            where each consecutive pair of real integrans is the real
            and imaginary parts of a single complex integrand, and you
            only care about the error in the complex plane rather than
            the error in the real and imaginary parts separately;
        - ERROR_L2
        - ERROR_L1
        - ERROR_LINF
            the absolute error is measured as |e| and the relative error
            as |err|/|val|, where |...| is the L1, L2, or L-infinity
            norm, respectively.  (|x| in the L1 norm is the sum of the
            absolute values of the components, in the L2 norm is the
            root mean square of the components, and in the L-infinity
            norm is the maximum absolute value of the components).
    maxEval : integer
        The maximum number of function evaluations.
    vectorized : boolean
        If ``vectorized=True`` the integration points are passed to the
        integrand function as an array of points, allowing parallel
        evaluation of different points.
    Returns
    -------
    val : numpy.ndarray
        The 1-D array of length ``fdim`` with the computed integral values
    err : numpy.ndarray
        The 1-D array of length ``fdim`` with the estimated errors. For
        smooth functions this estimate is usually conservative (see the
        results from the ``test_cubature.py`` script.
    Notes
    -----
    * The supplied function must return a 1-D ``np.ndarray`` object,
      even for single-valued functions
    References
    ----------
    .. [1] `Cubature (Multi-dimensional integration)
           <http://ab-initio.mit.edu/wiki/index.php/Cubature>`_.
    Examples
    --------
    >>> import numpy as np
    >>> from cubature import cubature
    Volume of a sphere:
    >>>
    >>> ndim = 3
    >>> xmin = np.zeros(ndim)
    >>> def integrand_sphere(x_array, *args):
    >>>     r, theta, phi = x_array
    >>>     return np.array([r**2*sin(phi)])
    >>>
    >>> radius = 1.
    >>> xmin = np.array([0, 0, 0], float)
    >>> xmax = np.array([radius, 2*pi, pi], float)
    >>> val, err = cubature(integrand_sphere, 3, xmin, xmax)
    More examples in ./examples/*.py
    """
    # checking xmin and xmax
    try:
        assert xmin.shape[0] == ndim
    except:
        raise ValueError('xmin.shape[0] is not equal ndim')
    try:
        assert xmax.shape[0] == ndim
    except:
        raise ValueError('xmax.shape[0] is not equal ndim')
    # getting fdim and checking function
    if vectorized:
        npt = 10
        dummy_x = np.ones((ndim*npt), dtype=np.float64)
        val = func(dummy_x, npt, *args)
        fdim = val.shape[0]//npt
    else:
        dummy_x = np.zeros((ndim), dtype=np.float64)
        val = func(dummy_x, *args)
        fdim = val.shape[0]
    try:
        assert val.ndim == 1
    except:
        raise ValueError('The input function must return a 1-D array!')

    method = _call_map.get((adaptive, vectorized), None)
    if method is None:
        s = 'unknown combination of adaptive (`{!r}`) and vectorized (`{!r}`).'
        s = s.format(adaptive, vectorized)
        raise ValueError(s)
    else:
        #val, err = _cython_cubature(func, ndim, xmin, xmax, method, abserr,
        #        relerr, norm, maxEval)
        pass
    
    #return val, err

#TODO
# - implement multiprocessing dividing the integration interval and spawning
# - a thread for each sub-interval...
# - perform a cProfile to see where the bottle nech actually is

