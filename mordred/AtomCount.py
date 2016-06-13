from rdkit.Chem import rdMolDescriptors

from ._base import Descriptor
from ._atomic_property import halogen

__all__ = (
    'AtomCount',
)


class AtomCount(Descriptor):
    r"""atom count descriptor.

    :type type: str
    :param type: type to count.

        * 'Atom'
        * 'HeavyAtom'
        * 'Spiro'
        * 'Bridgehead'
        * 'X' - all halogen
        * element symbol
    """
    __slots__ = ('_type',)

    @classmethod
    def preset(cls):
        return map(cls, [
            'Atom', 'HeavyAtom', 'Spiro', 'Bridgehead',
            'H', 'B', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'X',
        ])

    @property
    def explicit_hydrogens(self):
        u"""require explicit_hydrogens when type is 'H' or 'Atom'."""
        return self._type in set(['H', 'Atom'])

    def __str__(self):
        return 'n' + self._type

    def as_key(self):
        return self.__class__, (self._type,)

    def __init__(self, type='Atom'):
        self._type = type

    def _calc_X(self):
        X = halogen
        return sum(a.GetAtomicNum() in X for a in self.mol.GetAtoms())

    def _calc(self):
        return sum(a.GetSymbol() == self._type for a in self.mol.GetAtoms())

    def _calc_all(self):
        return self.mol.GetNumAtoms()

    def calculate(self):
        if self._type == 'X':
            return self._calc_X()
        elif self._type in ['Atom', 'HeavyAtom']:
            return self._calc_all()
        elif self._type == 'Spiro':
            return rdMolDescriptors.CalcNumSpiroAtoms(self.mol)
        elif self._type == 'Bridgehead':
            return rdMolDescriptors.CalcNumBridgeheadAtoms(self.mol)
        else:
            return self._calc()

    rtype = int
