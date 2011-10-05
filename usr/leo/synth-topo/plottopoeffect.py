import numpy
from matplotlib import pyplot
from fatiando import stats, gridder, logger, vis, potential
from fatiando.mesher.prism import Relief3D, relief2prisms, fill_relief


log = logger.get()
log.info(logger.header())

files = ['topo-ns', 'topo-ew', 'topo-diag', 'topo-mount']
#files = ['topo-ns.txt']
shape = (101, 101)
comps = ['gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']
for f in files:
    log.info("Plotting %s" % (f))
    pyplot.figure()
    pyplot.suptitle(f)
    for i, comp in enumerate(comps):
        xp, yp, zp, res = numpy.loadtxt('%s-%s.txt' % (f,comp),unpack=True,usecols=[0,1,2,3])
        pyplot.subplot(3,2,i+1)
        pyplot.title(comp)
        pyplot.axis('scaled')
        vis.pcolor(yp, xp, res, shape)
        pyplot.colorbar()
        vis.contour(yp, xp, -zp, shape, 5)
    pyplot.savefig('%s-ggt.png' % (f))
#pyplot.show()
