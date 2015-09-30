import numpy as np
from numpy import pi, sin

from cubature import cubature

def integrand_brick(x_array):
    return 1.

def integrand_sphere(x_array):
    r, theta, phi = x_array
    return r**2*sin(phi)

def integrand_ellipsoid(x_array, a, b, c):
    rho, phi, theta = x_array
    return a*b*c*rho**2*sin(theta)

def integrand_ellipsoid_v(x_array, a, b, c):
    rho = np.array(x_array[:, 0])
    phi = np.array(x_array[:, 1])
    theta = np.array(x_array[:, 2])
    return a*b*c*rho**2*sin(theta)

def exact_brick(a, b, c):
    return a*b*c

def exact_sphere(r):
    return 4./3*pi*r**3

def exact_ellipsoid(a,b,c):
    return 4./3*pi*a*b*c

if __name__ == '__main__':
    # brick
    print('_________________')
    print('')
    print('Brick')
    a, b, c = 1., 2., 3.
    xmin = np.zeros((3,), dtype=float)
    xmax = np.array([a, b, c], dtype=float)
    val, err = cubature(integrand_brick, 3, 1, xmin, xmax)
    print('Approximated: {0}'.format(val))
    print('Exact: {0}'.format(exact_brick(a,b,c)))
    # sphere
    print('_________________')
    print('')
    print('Sphere')
    radius = 1.
    xmin = np.array([0, 0, 0], np.float64)
    xmax = np.array([radius, 2*pi, pi], np.float64)
    val, err = cubature(integrand_sphere, 3, 1, xmin, xmax)
    print('Approximated: {0}'.format(val))
    print('Exact: {0}'.format(exact_sphere(radius)))
    # ellipsoid
    print('_________________')
    print('')
    print('Ellipsoid')
    a, b, c = 1., 2., 3.
    xmin = np.array([0, 0, 0], np.float64)
    xmax = np.array([1., 2*pi, pi], np.float64)
    val, err = cubature(integrand_ellipsoid, 3, 1, xmin, xmax, args=(a,b,c))
    print('Approximated: {0}'.format(val))
    print('Exact: {0}'.format(exact_ellipsoid(a,b,c)))
    print('_________________')
    # ellipsoid vectorized
    print('_________________')
    print('')
    print('Ellipsoid Vectorized')
    a, b, c = 1., 2., 3.
    xmin = np.array([0, 0, 0], np.float64)
    xmax = np.array([1., 2*pi, pi], np.float64)
    val, err = cubature(integrand_ellipsoid_v, 3, 1, xmin, xmax, args=(a,b,c),
        vectorized=True)
    print('Approximated: {0}'.format(val))
    print('Exact: {0}'.format(exact_ellipsoid(a,b,c)))
    print('_________________')

