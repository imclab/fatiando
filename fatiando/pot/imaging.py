"""
Imaging methods for potential fields.

Implements some of the methods described in Fedi and Pilkington (2012).
Most methods convert the observed data (gravity, magnetic, etc) into a physical
property distribution (density, magnetization, etc). Most methods require
gridded data to work.

* :func:`~fatiando.pot.imaging.geninv`: The Generalized Inverse solver in the
  frequency domain (Cribb, 1976)
* :func:`~fatiando.pot.imaging.sandwich`: Sandwich model (Pedersen, 1991). Uses
  depth weighting as in Pilkington (1997)
* :func:`~fatiando.pot.imaging.migrate`: 3D potential field migration
  (Zhdanov et al., 2011). Actually uses the formula of Fedi and Pilkington
  (2012), which are comprehensible.


**References**

Cribb, J. (1976), Application of the generalized linear inverse to the inversion
of static potential data, Geophysics, 41(6), 1365, doi:10.1190/1.1440686

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

Pedersen, L. B. (1991), Relations between potential fields and some equivalent
sources, Geophysics, 56(7), 961, doi:10.1190/1.1443129

Pilkington, M. (1997), 3-D magnetic imaging using conjugate gradients,
Geophysics, 62(4), 1132, doi:10.1190/1.1444214

Zhdanov, M. S., X. Liu, G. A. Wilson, and L. Wan (2011), Potential field
migration for rapid imaging of gravity gradiometry data, Geophysical
Prospecting, 59(6), 1052-1071, doi:10.1111/j.1365-2478.2011.01005.x

----
"""
import time

import numpy

from fatiando.msh.ddd import PrismMesh
from fatiando.pot import fourier
from fatiando.pot import prism as pot_prism
from fatiando.constants import G
from fatiando import utils
import fatiando.log

log = fatiando.log.dummy('fatiando.pot.imaging')


def migrate(x, y, z, gz, shape, zmin, zmax, nlayers, power=0.5, scale=1):
    """
    3D potential field migration (Zhdanov et al., 2011).

    Actually uses the formulation of Fedi and Pilkington (2012), which is
    comprehensible.

    .. warning:: Due to depth weighting, the values of the estimated physical
        property distribuition are out of scale. They should not be used as
        actual values of density or magnetization.

    .. note:: Only works on **gravity** data for now.

    .. note:: The data **must** be leveled, i.e., on the same height!

    .. note:: The coordinate system adopted is x->North, y->East, and z->Down

    Parameters:

    * x, y : 1D-arrays
        The x and y coordinates of the grid points
    * z : float or 1D-array
        The z coordinate of the grid points
    * gz : 1D-array
        The gravity anomaly data at the grid points
    * shape : tuple = (ny, nx)
        The shape of the grid
    * zmin, zmax : float
        The top and bottom, respectively, of the region where the physical
        property distribution is calculated
    * nlayers : int
        The number of layers used to divide the region where the physical
        property distribution is calculated
    * power : float
        The power law used for the depth weighting. This controls what depth
        the bulk of the solution will be.
    * scale : float
        A scale factor for the depth weights. Simply changes the scale of the
        physical property values.

    Returns:

    * mesh : :class:`fatiando.msh.ddd.PrismMesh`
        The estimated physical property distribution set in a prism mesh (for
        easy 3D plotting)

    """
    mesh = _makemesh(x, y, shape, zmin, zmax, nlayers)
    log.info("3D migration of gravity data:")
    if not isinstance(z, float) and not isinstance(z, int):
        z = z[0]
    else:
        z = z*numpy.ones_like(x)
    log.info("  data z coordinate: %g" % (z))
    log.info("  data shape: %s" % (str(shape)))
    log.info("  mesh zmin and zmax: %g, %g" % (zmin, zmax))
    log.info("  number of layers in the mesh: %d" % (nlayers))
    log.info("  depth weighting power law: %g" % (power))
    log.info("  depth weighting scale factor: %g" % (scale))
    tstart = time.clock()
    dx, dy, dz = mesh.dims
    # Synthetic tests show that its not good to offset the weights with the data
    # z coordinate. No idea why
    depths = mesh.get_zs()[:-1] + 0.5*dz
    weights = numpy.abs(depths)**power/(2*G*numpy.sqrt(numpy.pi))
    density = []
    for l in xrange(nlayers):
        sensibility_T = numpy.array(
            [pot_prism.gz(x, y, z, [p], dens=1) for p in mesh.get_layer(l)])
        density.extend(scale*weights[l]*numpy.dot(sensibility_T, gz))
    tend = time.clock()
    log.info("  total time for imaging: %s" % (utils.sec2hms(tend - tstart)))
    mesh.addprop('density', numpy.array(density))
    return mesh

def sandwich(x, y, z, data, shape, zmin, zmax, nlayers, power=0.5):
    """
    Sandwich model (Pedersen, 1991).

    Calculates a physical property distribution given potential field data on a
    **regular grid**. Uses depth weights.

    .. note:: Only works on **gravity** data for now.

    .. note:: The data **must** be leveled, i.e., on the same height!

    .. note:: The coordinate system adopted is x->North, y->East, and z->Down

    Parameters:

    * x, y : 1D-arrays
        The x and y coordinates of the grid points
    * z : float
        The z coordinate of the grid points
    * data : 1D-array
        The potential field at the grid points
    * shape : tuple = (ny, nx)
        The shape of the grid
    * zmin, zmax : float
        The top and bottom, respectively, of the region where the physical
        property distribution is calculated
    * nlayers : int
        The number of layers used to divide the region where the physical
        property distribution is calculated
    * power : float
        The power law used for the depth weighting. This controls what depth
        the bulk of the solution will be.

    Returns:

    * mesh : :class:`fatiando.msh.ddd.PrismMesh`
        The estimated physical property distribution set in a prism mesh (for
        easy 3D plotting)

    """
    mesh = _makemesh(x, y, shape, zmin, zmax, nlayers)
    log.info("Sandwich model of gravity data:")
    if not isinstance(z, float) and not isinstance(z, int):
        z = z[0]
    log.info("  data z coordinate: %g" % (z))
    log.info("  data shape: %s" % (str(shape)))
    log.info("  mesh zmin and zmax: %g, %g" % (zmin, zmax))
    log.info("  number of layers in the mesh: %d" % (nlayers))
    log.info("  depth weighting power law: %g" % (power))
    tstart = time.clock()
    freq, dataft = _getdataft(x, y, data, shape)
    dx, dy, dz = mesh.dims
    # Remove the last z because I only want depths to the top of the layers
    depths = mesh.get_zs()[:-1]
    weights = (numpy.abs(depths) + 0.5*dz)**(power)
    density = []
    # Offset by the data z because in the paper the data is at z=0
    for depth, weight in zip(depths - z, weights):
        density.extend(
            numpy.real(numpy.fft.ifft2(
                weight*(numpy.exp(-freq*depth) - numpy.exp(-freq*(depth + dz)))
                *freq*dataft /
                (numpy.pi*G*
                 reduce(numpy.add,
                     [w*(numpy.exp(-freq*h) - numpy.exp(-freq*(h + dz)))**2
                      + 10.**(-10) # To avoid zero division when freq[i]==0
                      for h, w in zip(depths, weights)])
                )
            ).ravel()))
    tend = time.clock()
    log.info("  total time for imaging: %s" % (utils.sec2hms(tend - tstart)))
    mesh.addprop('density', numpy.array(density))
    return mesh

def geninv(x, y, z, data, shape, zmin, zmax, nlayers):
    """
    Generalized Inverse imaging in the frequency domain (Cribb, 1976).

    Calculates a physical property distribution given potential field data on a
    **regular grid**.

    .. note:: Only works on **gravity** data for now.

    .. note:: The data **must** be leveled, i.e., on the same height!

    .. note:: The coordinate system adopted is x->North, y->East, and z->Down

    .. warning:: The Generalized Inverse does **not** use depth weights. This
        means that the solution will tend to be concentrated on the surface!

    Parameters:

    * x, y : 1D-arrays
        The x and y coordinates of the grid points
    * z : float
        The z coordinate of the grid points
    * data : 1D-array
        The potential field at the grid points
    * shape : tuple = (ny, nx)
        The shape of the grid
    * zmin, zmax : float
        The top and bottom, respectively, of the region where the physical
        property distribution is calculated
    * nlayers : int
        The number of layers used to divide the region where the physical
        property distribution is calculated

    Returns:

    * mesh : :class:`fatiando.msh.ddd.PrismMesh`
        The estimated physical property distribution set in a prism mesh (for
        easy 3D plotting)

    """
    mesh = _makemesh(x, y, shape, zmin, zmax, nlayers)
    log.info("Generalized Inverse imaging of gravity data:")
    if not isinstance(z, float) and not isinstance(z, int):
        z = z[0]
    log.info("  data z coordinate: %g" % (z))
    log.info("  data shape: %s" % (str(shape)))
    log.info("  mesh zmin and zmax: %g, %g" % (zmin, zmax))
    log.info("  number of layers in the mesh: %d" % (nlayers))
    tstart = time.clock()
    freq, dataft = _getdataft(x, y, data, shape)
    dx, dy, dz = mesh.dims
    # Remove the last z because I only want depths to the top of the layers
    depths = mesh.get_zs()[:-1] + 0.5*dz - z # Offset by the data height
    density = []
    for depth in depths:
        density.extend(
            numpy.real(
                numpy.fft.ifft2(
                    numpy.exp(-freq*depth)*freq*dataft/(numpy.pi*G)
                ).ravel()
            ))
    tend = time.clock()
    log.info("  total time for imaging: %s" % (utils.sec2hms(tend - tstart)))
    mesh.addprop('density', numpy.array(density))
    return mesh

def _getdataft(x, y, data, shape):
    """
    Get the Fourier transform of the data and the norm of the wavenumber vector
    """
    Fx, Fy = fourier._getfreqs(x, y, data, shape)
    freq = numpy.sqrt(Fx**2 + Fy**2)
    dataft = (2.*numpy.pi)*numpy.fft.fft2(numpy.reshape(data, shape))
    return freq, dataft

def _makemesh(x, y, shape, zmin, zmax, nlayers):
    """
    Make a prism mesh bounded by the data.
    """
    ny, nx = shape
    bounds = [x.min(), x.max(), y.min(), y.max(), zmin, zmax]
    mesh = PrismMesh(bounds, (nlayers, ny, nx))
    return mesh


