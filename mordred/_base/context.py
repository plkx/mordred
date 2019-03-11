from copy import copy

from rdkit import Chem

from .._util import conformer_to_numpy
from ..error import Missing3DCoordinate


class Context(object):
    __slots__ = "_mol_coords", "n_frags", "name", "_stack"

    def __init__(self, mol_coords, n_frags, name):
        self._mol_coords = mol_coords
        self.n_frags = n_frags
        self.name = name

    def __reduce_ex__(self, version):
        return self.__class__, (self._mol_coords, self.n_frags, self.name)

    def __str__(self):
        return self.name

    @classmethod
    def from_query(cls, mol, mol_states, id):
        if not isinstance(mol, Chem.Mol):
            raise TypeError("{!r} is not rdkit.Chem.Mol instance".format(mol))

        n_frags = len(Chem.GetMolFrags(mol))

        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        else:
            name = Chem.MolToSmiles(Chem.RemoveHs(mol, updateExplicitCount=True))

        mol_coords = {}
        num_imp_hs = sum(a.GetNumImplicitHs() for a in mol.GetAtoms())

        for state in mol_states:
            eh, ke, r3d = state
            added = False
            if eh and num_imp_hs > 0:
                added = True
                m = Chem.AddHs(mol)
            elif not eh:
                m = Chem.RemoveHs(mol, updateExplicitCount=True)
            else:
                m = copy(mol)

            coord = None
            try:
                conf = m.GetConformer(id)
            except ValueError:
                pass
            else:
                if r3d and not added and conf.Is3D():
                    coord = conformer_to_numpy(conf)

            if ke:
                Chem.Kekulize(m)
            else:
                Chem.SanitizeMol(m)

            m.RemoveAllConformers()
            mol_coords[state] = m, coord

        return cls(mol_coords, n_frags, name)

    @classmethod
    def from_calculator(cls, calc, mol, id):
        return cls.from_query(mol, calc._mol_states, id)

    def get_coord(self, desc):
        _, crd = self._mol_coords[desc._get_state_key()]
        if crd is None:
            return desc.fail(Missing3DCoordinate())

        return crd

    def get_mol(self, desc):
        mol, _ = self._mol_coords[desc._get_state_key()]
        return mol

    def reset(self):
        self._stack = []

    def add_stack(self, d):
        self._stack.append(d)

    def get_stack(self):
        return self._stack
