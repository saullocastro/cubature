import pstats, cProfile

import pyximport
pyximport.install()

import numpy as np
from cubature import cubature
import cubature._test_integrands as ti

u = 2*np.pi*15/609
a = np.array([19.51, 2], dtype=float)
n = 1

xmin = np.zeros((2,))
xmax = np.ones((2,)) * n

val, err = cubature(ti.genz_oscillatory, 2, 1, xmin, xmax, args=(a, u),
        adaptive='p')

cProfile.runctx("cubature(ti.genz_oscillatory, 2, 1, xmin, xmax, args=(a,u), adaptive='h')", globals(), locals(), "cubature.prof")

s = pstats.Stats("cubature.prof")
s.strip_dirs().sort_stats("time").print_stats()
