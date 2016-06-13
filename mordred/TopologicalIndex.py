from __future__ import division

from ._base import Descriptor
from ._graph_matrix import Radius as CRadius
from ._graph_matrix import Diameter as CDiameter

__all__ = ('Diameter', 'Radius', 'TopologicalShapeIndex', 'PetitjeanIndex',)


class TopologicalIndexBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    def as_key(self):
        return self.__class__, ()

    rtype = int


class Radius(TopologicalIndexBase):
    r"""radius descriptor."""
    __slots__ = ()

    def __str__(self):
        return 'Radius'

    def dependencies(self):
        return {'R': CRadius(self.explicit_hydrogens)}

    def calculate(self, R):
        return int(R)


class Diameter(TopologicalIndexBase):
    r"""diameter descriptor."""
    __slots__ = ()

    def __str__(self):
        return 'Diameter'

    def dependencies(self):
        return {'D': CDiameter(self.explicit_hydrogens)}

    def calculate(self, D):
        return int(D)


class TopologicalShapeIndex(TopologicalIndexBase):
    r"""topological shape index descriptor.

    .. math::

        I_{\rm topo} = \frac{D - R}{R}

    where
    :math:`R` is graph radius,
    :math:`D` is graph diameter.

    :returns: NaN when :math:`R = 0`
    """
    __slots__ = ()

    def __str__(self):
        return 'TopoShapeIndex'

    def dependencies(self):
        return {
            'R': CRadius(self.explicit_hydrogens),
            'D': CDiameter(self.explicit_hydrogens),
        }

    def calculate(self, R, D):
        with self.rethrow_zerodiv():
            return (D - R) / R

    rtype = float


class PetitjeanIndex(TopologicalShapeIndex):
    r"""Petitjean index descriptor.

    .. math::

        I_{\rm Petitjean} = \frac{D - R}{D}

    where
    :math:`R` is graph radius,
    :math:`D` is graph diameter.

    :returns: NaN when :math:`D = 0`
    """
    __slots__ = ()

    def __str__(self):
        return 'PetitjeanIndex'

    def calculate(self, R, D):
        with self.rethrow_zerodiv():
            return (D - R) / D
