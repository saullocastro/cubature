import pstats, cProfile

import pyximport
pyximport.install()

import numpy as np
from cubature import cubature
import cubature._test_integrands as ti

u = 2*np.pi*15/609
a = np.array([19.51, 2], dtype=float)
n = 10 

xmin = np.zeros((2,))
xmax = np.ones((2,)) * n

cProfile.runctx("cubature(ti.genz_oscillatory, 2, 1, xmin, xmax, abserr=1e-4, relerr=1e-8, args=(a,u), adaptive='h')", globals(), locals(), "cubature_1.prof")

s = pstats.Stats("cubature_1.prof")
s.strip_dirs().sort_stats("time").print_stats()

cProfile.runctx("cubature(ti.genz_oscillatory_c, 2, 1, xmin, xmax, abserr=1e-4, relerr=1e-8, args=(a,u), adaptive='h')", globals(), locals(), "cubature_2.prof") 

s = pstats.Stats("cubature_2.prof")
s.strip_dirs().sort_stats("time").print_stats()

