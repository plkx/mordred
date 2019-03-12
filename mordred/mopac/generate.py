from rdkit import Chem
import numpy as np


def distance(i, j):
    return np.linalg.norm(i - j)


def angle(i, j, k):
    rij = i - j
    rkj = k - j
    return np.arccos(np.clip(rij.dot(rkj) / (np.linalg.norm(rij) * np.linalg.norm(rkj)), -1, 1)) / np.pi * 180


def dihedral(i, j, k, l):
    rij = i - j
    rkj = k - j
    rkl = k - l
    cijk = np.cross(rij, rkj)
    cjkl = np.cross(rkj, rkl)
    s = np.sign(rkj.dot(np.cross(cijk, cjkl)))
    return s * np.arccos(np.clip(cijk.dot(cjkl) / (np.linalg.norm(cijk) * np.linalg.norm(cjkl)), -1, 1)) / np.pi * 180


def internal(mol, confId=-1):
    conf = mol.GetConformer(confId)
    it = zip(mol.GetAtoms(), (np.array(conf.GetAtomPosition(i)) for i in range(conf.GetNumAtoms())))
    atm1, crd1 = next(it)
    yield atm1, (0., 0., 0.), (0, 0, 0)
    atm2, crd2 = next(it)
    yield atm2, (distance(crd2, crd1), 0., 0.), (1, 0, 0)
    atm3, crd3 = next(it)
    yield atm3, (distance(crd3, crd1), angle(crd3, crd1, crd2), 0.), (1, 2, 0)

    for atmN, crdN in it:
        yield atmN, (distance(crdN, crd1), angle(crdN, crd1, crd2), dihedral(crdN, crd1, crd2, crd3)), (1, 2, 3)


def generate_mopac_cartesian_input(mol, fp, condition="PM3 MMOK", confId=-1):
    fc = Chem.GetFormalCharge(mol)
    conf = mol.GetConformer(confId)
    with fp:
        fp.write(" ".join(condition.split("\n")) + " XYZ CHARGE={}\n\n\n".format(fc))
        for i in range(mol.GetNumAtoms()):
            a = mol.GetAtomWithIdx(i)
            f = int(not (a.GetBoolProp("Fixed") if a.HasProp("Fixed") else False))
            x, y, z = conf.GetAtomPosition(i)
            fp.write("{} {} {} {} {} {} {}\n".format(a.GetSymbol(), x, f, y, f, z, f))


def generate_mopac_internal_input(mol, fp, condition="PM3 MMOK", confId=-1):
    fc = Chem.GetFormalCharge(mol)
    with fp:
        fp.write(" ".join(condition.split("\n")) + " XYZ CHARGE={}\n\n\n".format(fc))
        for a, (dist, ang, dihe), (i, j, k) in internal(mol, confId):
            f = int(not (a.GetBoolProp("Fixed") if a.HasProp("Fixed") else False))
            fp.write("{} {:20.14f} {} {:20.14f} {} {:20.14f} {} {} {} {}\n".format(
                a.GetSymbol(),
                dist, f, ang, f, dihe, f,
                i, j, k,
            ))


def main():
    from rdkit.Chem import AllChem as Chem
    import sys

    mol = Chem.AddHs(Chem.MolFromSmiles(sys.stdin.read()))
    Chem.EmbedMolecule(mol)
    Chem.MMFFOptimizeMolecule(mol)
    generate_mopac_cartesian_input(mol, sys.stdout)
    generate_mopac_internal_input(mol, sys.stderr)


if __name__ == "__main__":
    main()
