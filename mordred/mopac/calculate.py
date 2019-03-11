import os
from tempfile import TemporaryDirectory
import subprocess

from .generate import generate_mopac_input
from .parser import parse_output, get_dipole_from_arc


def get_mopac_path():
    return 'run_mopac7'


def calculate(mol, confId=-1, condition='PM3 MMOK', executable=None):
    if executable is None:
        executable = get_mopac_path()

    with TemporaryDirectory() as d:
        with open(os.path.join(d, "mol.dat"), 'w') as i:
            generate_mopac_input(mol, i, condition=condition, confId=confId)
            subprocess.run([executable, 'mol'], cwd=d, stdout=subprocess.PIPE)
            result = parse_output(open(os.path.join(d, "mol.OUT"), 'r'))
            result['dipole'] = get_dipole_from_arc(open(os.path.join(d, "mol.arc"), 'r'))

    conf = mol.GetConformer(confId)
    for i, xyz in enumerate(result["coordinates"]):
        conf.SetAtomPosition(i, xyz)

    return result


def main():
    import sys
    from rdkit.Chem import AllChem as Chem
    from contextlib import closing

    mol = Chem.AddHs(Chem.MolFromSmiles(sys.stdin.read()))
    Chem.EmbedMolecule(mol)
    Chem.MMFFOptimizeMolecule(mol)

    with closing(Chem.SDWriter(sys.stdout)) as w:
        w.write(mol)

        result = calculate(mol)
        print(result['heat_of_formation'], file=sys.stderr)
        result = calculate(mol, condition="PM3 MMOK PRECISE")
        print(result['heat_of_formation'], file=sys.stderr)
        mol.SetDoubleProp("heat of formation", result['heat_of_formation'])
        mol.SetDoubleProp("total energy", result['total_energy'])
        mol.SetDoubleProp("electronic energy", result['electronic_energy'])
        mol.SetDoubleProp("core-core repulsion", result['core_core_repulsion'])
        mol.SetDoubleProp("ionization potential", result['ionization_potential'])
        mol.SetDoubleProp("dipole", result['dipole'])
        fl = result['num_of_filled_levels']
        homo = max(result['eigenvalues'][:fl])
        lumo = min(result['eigenvalues'][fl:])
        mol.SetDoubleProp("HOMO", homo)
        mol.SetDoubleProp("LUMO", lumo)
        mol.SetDoubleProp("HOMO-LUMO gap", lumo - homo)

        w.write(mol)


if __name__ == "__main__":
    main()
