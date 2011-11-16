"""
Generate synthetic topographies and save them to files
"""
from math import sin, cos, pi
import numpy
from matplotlib import pyplot
from fatiando import stats, gridder, logger, vis

log = logger.get()
log.info(logger.header())

log.info("Generating synthetic topography")
area = (-7000, 7000, -7000, 7000)
shape = (141, 141)
x, y = gridder.regular(area, shape)

topons = (200 + 200*stats.gaussian2d(x, y, 0, 0, cov=[[50000**2,0],[0,1000**2]]))
numpy.savetxt('topo-ns.txt', numpy.array([x, y, topons]).T)

topoew = (200 + 200*stats.gaussian2d(x, y, 0, 0, cov=[[1000**2,0],[0,50000**2]]))
numpy.savetxt('topo-ew.txt', numpy.array([x, y, topoew]).T)

sigx, sigy = 50000.**2, 1000.**2
theta = -45.*pi/180.
a = sigx*cos(theta)**2 + sigy*sin(theta)**2
b = (sigy - sigx)*cos(theta)*sin(theta)
c = sigx*sin(theta)**2 + sigy*cos(theta)**2
topodiag = (200 + 200*stats.gaussian2d(x, y, 0, 0, cov=[[a,b],[b,c]]))
numpy.savetxt('topo-diag.txt', numpy.array([x, y, topodiag]).T)

topomount = (200 + 200*stats.gaussian2d(x, y, 0, 0, cov=[[2000**2,0],[0,2000**2]]))
numpy.savetxt('topo-mount.txt', numpy.array([x, y, topomount]).T)

pyplot.figure()
pyplot.suptitle("Synthetic topography")
pyplot.subplot(2,2,1)
pyplot.title("NS")
pyplot.axis('scaled')
vis.pcolor(y*0.001, x*0.001, topons, shape)
pyplot.colorbar()
pyplot.subplot(2,2,2)
pyplot.title("EW")
pyplot.axis('scaled')
vis.pcolor(y*0.001, x*0.001, topoew, shape)
pyplot.colorbar()
pyplot.subplot(2,2,3)
pyplot.title("Diagonal")
pyplot.axis('scaled')
vis.pcolor(y*0.001, x*0.001, topodiag, shape)
pyplot.colorbar()
pyplot.subplot(2,2,4)
pyplot.title("Mount")
pyplot.axis('scaled')
vis.pcolor(y*0.001, x*0.001, topomount, shape)
pyplot.colorbar()
pyplot.savefig('topo.png')
pyplot.show()


