"""
Generate synthetic topographies and save them to files
"""
import numpy
from matplotlib import pyplot
from fatiando import stats, gridder, logger, vis, potential
from fatiando.mesher.prism import Relief3D, relief2prisms, fill_relief


log = logger.get()
log.info(logger.header())

#files = ['topo-ns.txt', 'topo-ew.txt', 'topo-diag.txt', 'topo-mount.txt']
files = ['topo-ns.txt']
shape = (101, 101)
area = (-5000, 5000, -5000, 5000)
areatopo = (-7000, 7000, -7000, 7000)
shapetopo = (141, 141)
dens = 1000.
offset = 150.
comps = {'gxx':potential.prism.gxx, 'gxy':potential.prism.gxy,
         'gxz':potential.prism.gxz, 'gyy':potential.prism.gyy,
         'gyz':potential.prism.gyz, 'gzz':potential.prism.gzz}

for f in files:
    log.info("Load topography from %s" % (f))
    x, y, height = numpy.loadtxt(f, unpack=True)
    scalars = [dens]*len(height)
    nodes = (x,y,-1*height)
    relief = fill_relief(scalars,Relief3D(0,gridder.spacing(areatopo,shapetopo),nodes))

    xp, yp, tmp = gridder.cut(x, y, [height], -5000., 5000., -5000., 5000.)
    zp = -1.*tmp[0] - offset
    log.info("Calculate tensor")
    pyplot.figure()
    pyplot.suptitle(f)
    for i, comp in enumerate(comps):
        res = comps[comp](xp, yp, zp, relief2prisms(relief, 'density'))
        numpy.savetxt('%s-%s.txt' % (f, comp), numpy.array([xp,yp,zp,res]).T)
        pyplot.subplot(3,2,i+1)
        pyplot.title(comp)
        pyplot.axis('scaled')
        vis.pcolor(yp, xp, res, shape)
        pyplot.colorbar()
pyplot.show()
