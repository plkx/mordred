from __future__ import division

import numpy as np
from rdkit import Chem
from networkx import Graph, floyd_warshall_numpy

from ._base import Descriptor
from ._atomic_property import AtomicProperty, get_properties
from ._matrix_attributes import methods, get_method

__all__ = ('BaryszMatrix',)


class BaryszMatrixBase(Descriptor):
    explicit_hydrogens = False
    __slots__ = ()


class Barysz(BaryszMatrixBase):
    __slots__ = ('_prop',)

    _carbon = Chem.Atom(6)

    def as_key(self):
        return self.__class__, (self._prop,)

    def __init__(self, prop):
        self._prop = prop

    def dependencies(self):
        return {'P': self._prop}

    def calculate(self, P):
        C = self._prop.prop(self._carbon)

        G = Graph()

        G.add_nodes_from(a.GetIdx() for a in self.mol.GetAtoms())

        for bond in self.mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            pi = bond.GetBondTypeAsDouble()

            w = (C * C) / (P[i] * P[j] * pi)

            G.add_edge(i, j, weight=w)

        sp = floyd_warshall_numpy(G)
        np.fill_diagonal(sp, [1. - C / P[a.GetIdx()] for a in self.mol.GetAtoms()])
        return sp


class BaryszMatrix(BaryszMatrixBase):
    r"""barysz matrix descriptor.

    :type prop: :py:class:`str` or :py:class:`function`
    :param prop: :ref:`atomic_properties`

    :type type: str
    :param type: :ref:`matrix_aggregating_methods`

    :returns: NaN when any properties are NaN
    """
    __slots__ = ('_prop', '_type',)

    @classmethod
    def preset(cls):
        return (cls(p, m) for p in get_properties() for m in methods)

    def __str__(self):
        return '{}_Dz{}'.format(self._type.__name__, self._prop.as_argument)

    def as_key(self):
        return self.__class__, (self._prop, self._type)

    def __init__(self, prop='Z', type='SpMax'):
        self._prop = AtomicProperty(self.explicit_hydrogens, prop)
        self._type = get_method(type)

    def dependencies(self):
        return dict(
            result=self._type(
                Barysz(self._prop),
                self.explicit_hydrogens,
                self.kekulize,
            )
        )

    def calculate(self, result):
        return result

    rtype = float
