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
from numpy import log, arctan2, sqrt

SI2EOTVOS = 1000000000.0
SI2MGAL = 100000.0
G = 0.00000000006673

def gz(xp, yp, zp, prisms):
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
    res = numpy.zeros(len(xp))
    for prism in prisms:
        density = prism.props['density']
        # First thing to do is make the computation point P the origin of the
        # coordinate system
        x = [prism.x1 - xp, prism.x2 - xp]
        y = [prism.y1 - yp, prism.y2 - yp]
        z = [prism.z1 - zp, prism.z2 - zp]
        # Evaluate the integration limits 
        for k in range(2):
            for j in range(2):
                for i in range(2):
                    r = sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
                    kernel = (x[i]*log(y[j] + r) +
                              y[j]*log(x[i] + r) -
                              z[k]*arctan2(x[i]*y[j], z[k]*r))
                    res += ((-1)**(i + j + k))*kernel*density
    # Now all that is left is to multiply res by the gravitational constant and
    # convert it to mGal units 
    res *= G*SI2MGAL;
    return res
