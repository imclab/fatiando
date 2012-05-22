#import pyximport; pyximport.install()
#import prism
#from fatiando.mesher.ddd import Prism
#from fatiando import gridder, vis, potential
import timeit

setup = """
import pyximport; pyximport.install()
import cprism
import prism
from fatiando.mesher.ddd import Prism
from fatiando import gridder, vis, potential
prisms = [Prism(-4000,-3000,-4000,-3000,0,2000,{'density':1000}),
          Prism(-1000,1000,-1000,1000,0,2000,{'density':-1000}),
          Prism(2000,4000,3000,4000,0,2000,{'density':200})]
shape = (500,500)
xp, yp, zp = gridder.regular((-5000, 5000, -5000, 5000), shape, z=-100)
"""
n = 20
print "potential.prism:", timeit.timeit("potential.prism.gz(xp, yp, zp, prisms)", setup, number=n)/float(n)
print "numpy prism:", timeit.timeit("prism.gz(xp, yp, zp, prisms)", setup, number=n)/float(n)
print "cython prism:", timeit.timeit("cprism.gz(xp, yp, zp, prisms)", setup, number=n)/float(n)
#
#vis.map.pyplot.figure()
#vis.map.pyplot.axis('scaled')
#vis.map.contour(xp, yp, newgz, shape, 12, color='r')
#vis.map.contour(xp, yp, gz, shape, 12, color='k', style='dashed')
#vis.map.pyplot.show()
