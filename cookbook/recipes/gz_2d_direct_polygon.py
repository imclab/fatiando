"""
2D direct gravity modeling with polygons
"""
from matplotlib import pyplot
import numpy
from fatiando import potential, vis, logger
from fatiando.mesher.ddd import draw_polygon
from fatiando.mesher.dd import Polygon

log = logger.get()
log.info(logger.header())
log.info(__doc__)

# Notice that the last two number are switched.
# This way, the z axis in the plots points down.
area = (-5000, 5000, 5000, 0)
axes = pyplot.figure().add_subplot(1,1,1)
pyplot.xlabel("X")
pyplot.ylabel("Z")
polygons = [Polygon(draw_polygon(area, axes), {'density':500})]
xp = numpy.arange(-4500, 4500, 100)
zp = numpy.zeros_like(xp)
gz = potential.talwani.gz(xp, zp, polygons)

pyplot.figure()
pyplot.axis('scaled')
pyplot.subplot(2,1,1)
pyplot.title(r"Gravity anomaly produced by the model")
pyplot.plot(xp, gz, '.-k')
pyplot.ylabel("mGal")
pyplot.xlim(-5000, 5000)
pyplot.subplot(2,1,2)
vis.map.polygon(polygons[0], '.-k')
pyplot.xlabel("X")
pyplot.ylabel("Z")
vis.map.set_area(area)
pyplot.show()
