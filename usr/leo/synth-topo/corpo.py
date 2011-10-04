
import numpy
from matplotlib import pyplot
from fatiando import stats, gridder, logger, vis, potential
from fatiando.mesher.prism import Prism3D


log = logger.get()
log.info(logger.header())

prism = Prism3D(-1500,1500,-1500,1500,-200,0,{'density':1100.})

shape = (101, 101)

comps = {'gxx':potential.prism.gxx, 'gxy':potential.prism.gxy,
         'gxz':potential.prism.gxz, 'gyy':potential.prism.gyy,
         'gyz':potential.prism.gyz, 'gzz':potential.prism.gzz}

files = ['topo-ns', 'topo-ew', 'topo-diag', 'topo-mount']
for f in files:
    log.info("Calculate tensor for %s" % (f))
    xp, yp, zp = numpy.loadtxt('%s-gxx.txt' % (f), unpack=True, usecols=[0,1,2])
    pyplot.figure()
    pyplot.suptitle('corpo %s' % (f))
    for i, comp in enumerate(comps):
        res = comps[comp](xp, yp, zp, [prism])
        numpy.savetxt('%s-corpo-%s.txt' % (f, comp), numpy.array([xp,yp,zp,res]).T)
        pyplot.subplot(3,2,i+1)
        pyplot.title(comp)
        pyplot.axis('scaled')
        vis.pcolor(yp, xp, res, shape)
        pyplot.colorbar()
    pyplot.savefig('%s-corpo.png' % (f))
