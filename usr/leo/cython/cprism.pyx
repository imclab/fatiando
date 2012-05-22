"""
Calculate the potential fields and derivatives of the 3D right rectangular prism 
using the forumla of Nagy et al. (2000)

**Gravity**

* :func:`~fatiando.potential.prism.gz`
* :func:`~fatiando.potential.prism.gxx`
* :func:`~fatiando.potential.prism.gxy`
* :func:`~fatiando.potential.prism.gxz`
* :func:`~fatiando.potential.prism.gyy`
* :func:`~fatiando.potential.prism.gyz`
* :func:`~fatiando.potential.prism.gzz`

**Magnetic**


**References**

Nagy, D., G. Papp, and J. Benedek, 2000, The gravitational potential and its
derivatives for the prism: Journal of Geodesy, 74, 552--560,
doi: 10.1007/s001900000116.
    
:author: Leonardo Uieda (leouieda@gmail.com)
:date: Created 11-Sep-2010
:license: GNU Lesser General Public License v3 (http://www.gnu.org/licenses/)

----

"""
import numpy
from libc.math cimport log, atan2, sqrt
# Import Cython definitions for numpy
cimport numpy

DTYPE = numpy.float
ctypedef numpy.float_t DTYPE_T

SI2EOTVOS = 1000000000.0
SI2MGAL = 100000.0
G = 0.00000000006673

def gz(numpy.ndarray[DTYPE_T, ndim=1] xp, numpy.ndarray[DTYPE_T, ndim=1] yp,
       numpy.ndarray[DTYPE_T, ndim=1] zp, prisms):
    """
    Calculates the :math:`g_z` gravity acceleration component.

    .. note:: The coordinate system of the input parameters is to be x -> North,
        y -> East and z -> **DOWN**.

    .. note:: All input values in **SI** units(!) and output in **mGal**!

    Parameters:
    
    * xp, yp, zp : arrays
        Arrays with the x, y, and z coordinates of the computation points.
    * prisms : list of :class:`~fatiando.mesher.ddd.Prism`
        The density model used to calculate the gravitational effect.
        Prisms must have the property ``'density'``. Prisms that don't have this
        property will be ignored in the computations. Elements of *prisms* that
        are None will also be ignored. *prisms* can also be a
        :class:`~fatiando.mesher.ddd.PrismMesh`.

    Returns:
    
    * gz : array
        The :math:`g_z` component calculated on *points*

    """
    cdef int i, j, k, l, size
    cdef numpy.ndarray[DTYPE_T, ndim=1] res
    cdef DTYPE_T density, kernel, r
    cdef DTYPE_T x1, x2, y1, y2, z1, z2, dx1, dx2, dy1, dy2, dz1, dz2
    size = len(xp)
    res = numpy.zeros(size, dtype=DTYPE)
    for prism in prisms:
        density = prism.props['density']
        x1, x2 = prism.x1, prism.x2
        y1, y2 = prism.y1, prism.y2
        z1, z2 = prism.z1, prism.z2
        for l in xrange(size):
            # First thing to do is make the computation point P the origin of the
            # coordinate system
            dx1, dx2 = x1 - xp[l], x2 - xp[l]
            dy1, dy2 = y1 - yp[l], y2 - yp[l]
            dz1, dz2 = z1 - zp[l], z2 - zp[l]
            # Evaluate the integration limits 
            r = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
            kernel += (dx1*log(dy1 + r) + dy1*log(dx1 + r) -
                       dz1*atan2(dx1*dy1, dz1*r))
            r = sqrt(dx2*dx2 + dy1*dy1 + dz1*dz1)
            kernel += -(dx2*log(dy1 + r) + dy1*log(dx2 + r) -
                        dz1*atan2(dx2*dy1, dz1*r))
            r = sqrt(dx1*dx1 + dy2*dy2 + dz1*dz1)
            kernel += -(dx1*log(dy2 + r) + dy2*log(dx1 + r) -
                        dz1*atan2(dx1*dy2, dz1*r))
            r = sqrt(dx2*dx2 + dy2*dy2 + dz1*dz1)
            kernel += (dx2*log(dy2 + r) + dy2*log(dx2 + r) -
                       dz1*atan2(dx2*dy2, dz1*r))
            r = sqrt(dx1*dx1 + dy1*dy1 + dz2*dz2)
            kernel += -(dx1*log(dy1 + r) + dy1*log(dx1 + r) -
                        dz2*atan2(dx1*dy1, dz2*r))
            r = sqrt(dx2*dx2 + dy1*dy1 + dz2*dz2)
            kernel += (dx2*log(dy1 + r) + dy1*log(dx2 + r) -
                       dz2*atan2(dx2*dy1, dz2*r))
            r = sqrt(dx1*dx1 + dy2*dy2 + dz2*dz2)
            kernel += (dx1*log(dy2 + r) + dy2*log(dx1 + r) -
                       dz2*atan2(dx1*dy2, dz2*r))
            r = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)
            kernel += -(dx2*log(dy2 + r) + dy2*log(dx2 + r) -
                        dz2*atan2(dx2*dy2, dz2*r))
            res[l] += kernel*density
    # Now all that is left is to multiply res by the gravitational constant and
    # convert it to mGal units 
    res *= G*SI2MGAL
    return res
