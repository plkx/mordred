from .mopac.calculate import calc_energy
from ._base.descriptor import Descriptor

__all__ = ("Mopac",)


class MopacCache(Descriptor):
    since = "1.2.0"
    require_3D = True
    explicit_hydrogens = True

    __slots__ = ("_hamiltonian", "_timeout")

    def __init__(self, hamiltonian=None, timeout=None):
        self._hamiltonian = hamiltonian
        self._timeout = timeout

    def parameters(self):
        return (self._hamiltonian, self._timeout)

    def calculate(self):
        return calc_energy(self.get_3D_mol(), hamiltonian=self._hamiltonian, timeout=self._timeout)


type_to_attr = {
    "dipole": "dipole",
    "E": "total_energy",
    "Eele": "electronic_energy",
    "HF": "heat_of_formation",
    "IP": "ionization_potential",
    "LUMO": "lumo",
    "HOMO": "homo",
}


class Mopac(Descriptor):
    types = {
        "dipole": "dipole moment",
        "E": "total energy (kcal/mol)",
        "Eele": "electronic energy (kcal/mol)",
        "HF": "heat of formation (kcal/mol)",
        "IP": "ionization potential (kcal/mol)",
        "LUMO": "energy (eV) of the Lowest Unoccupied Molecular Orbital",
        "HOMO": "energy (eV) of the Highest Occupied Molecular Orbital",
    }

    since = "1.2.0"
    require_3D = True
    explicit_hydrogens = True

    __slots__ = ("_hamiltonian", "_type", "_timeout")

    @classmethod
    def preset(cls, version):
        return (cls(h, t) for h in ["AM1", "MNDO", "PM3"] for t in cls.types)

    def description(self):
        return "The {} calculated usin the {} hamiltonian".format(
            self.types[self._type], self._hamiltonian)

    def __str__(self):
        return "{}_{}".format(self._hamiltonian, self._type)

    def parameters(self):
        return (self._hamiltonian, self._type, self._timeout)

    def dependencies(self):
        return {"result": MopacCache(self._hamiltonian, self._timeout)}

    def calculate(self, result):
        return getattr(result, type_to_attr[self._type])

    def __init__(self, hamiltonian="AM1", type="total_energy", timeout=60):
        self._hamiltonian = hamiltonian
        self._type = type
        self._timeout = timeout

    rtype = float
