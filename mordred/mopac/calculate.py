from __future__ import print_function

import os
import sys
import time
import subprocess
from tempfile import TemporaryDirectory

from rdkit.Chem import rdForceFieldHelpers

from .parser import parse_output, get_dipole_from_arc
from .generate import generate_mopac_input


def get_mopac_path():
    if sys.platform.startswith("win"):
        return "MOPAC_7.1.exe"
    else:
        return "run_mopac7"


def calculate(mol, condition="PM3 XYZ MMOK", confId=-1, timeout=None, executable=None):
    if sum(a.GetNumImplicitHs() for a in mol.GetAtoms()) != 0:
        raise ValueError("mol has implicit hydrogens")

    if executable is None:
        executable = get_mopac_path()

    with TemporaryDirectory() as d:
        with open(os.path.join(d, "mol.dat"), "w") as i:
            generate_mopac_input(mol, i, condition=condition, confId=confId)
            subprocess.run(
                [executable, "mol"], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            result = parse_output(open(os.path.join(d, "mol.OUT"), "r"))
            result.dipole = get_dipole_from_arc(open(os.path.join(d, "mol.arc"), "r"))

    conf = mol.GetConformer(confId)
    for i, xyz in enumerate(result.coordinates):
        conf.SetAtomPosition(i, xyz)

    return result


def calc_energy(mol, hamiltonian="PM3", confId=-1, timeout=None, executable=None):
    return calculate(mol, "{} XYZ MMOK 1SCF".format(hamiltonian), confId, timeout, executable)


def optimize(mol, hamiltonian="PM3", confId=-1, timeout=None, executable=None):
    t1 = time.time()
    rdForceFieldHelpers.MMFFOptimizeMolecule(mol, maxIters=10000, confId=confId)
    t2 = time.time()
    calculate(
        mol,
        "{} XYZ MMOK".format(hamiltonian),
        confId,
        timeout - (t2 - t1) if timeout is not None else None,
        executable,
    )
    t3 = time.time()
    return calculate(
        mol,
        "{} XYZ MMOK PRECISE".format(hamiltonian),
        confId,
        timeout - (t3 - t1) if timeout is not None else None,
        executable,
    )


def main():
    import sys
    from rdkit.Chem import AllChem as Chem
    from contextlib import closing

    mol = Chem.AddHs(Chem.MolFromSmiles(sys.stdin.read()))
    Chem.EmbedMolecule(mol)

    with closing(Chem.SDWriter(sys.stdout)) as w:
        w.write(mol)

        result = optimize(mol)
        print(result.heat_of_formation, file=sys.stderr)  # noqa: T001
        mol.SetDoubleProp("heat of formation", result.heat_of_formation)
        mol.SetDoubleProp("total energy", result.total_energy)
        mol.SetDoubleProp("electronic energy", result.electronic_energy)
        mol.SetDoubleProp("core-core repulsion", result.core_core_repulsion)
        mol.SetDoubleProp("ionization potential", result.ionization_potential)
        mol.SetDoubleProp("dipole", result.dipole)
        mol.SetDoubleProp("HOMO", result.homo)
        mol.SetDoubleProp("LUMO", result.lumo)
        mol.SetDoubleProp("HOMO-LUMO gap", result.lumo - result.homo)

        w.write(mol)


if __name__ == "__main__":
    main()
