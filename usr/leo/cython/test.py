import timeit

setup = """
import pyximport; pyximport.install()
import cprism
import prism
import fatiando as ft
prisms = ft.msh.ddd.PrismMesh((-2000,2000,-2000,2000,0,1000), (10,10,10))
prisms.addprop('density', [1000]*prisms.size)
shape = (100, 100)
xp, yp, zp = ft.grd.regular((-5000, 5000, -5000, 5000), shape, z=-100)
"""
n = 3
print "potential.prism:", timeit.timeit("ft.pot.prism.gz(xp, yp, zp, prisms)", setup, number=n)/float(n)
print "numpy prism:", timeit.timeit("prism.gz(xp, yp, zp, prisms)", setup, number=n)/float(n)
print "cython prism:", timeit.timeit("cprism.gz(xp, yp, zp, prisms)", setup, number=n)/float(n)
