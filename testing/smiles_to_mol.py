from rdkit import Chem
from rdkit.Chem import AllChem

SMILES = "C1(CCCCCC2CCCCCC3CCC4CCCCCC(CCCCC5)CC2)CCCCCC(CCCCC3)CCC(CCCCC4)CCCCCC5CC1"

mol_noh = Chem.MolFromSmiles(SMILES)
mol_withh = Chem.AddHs(mol_noh)
Chem.AllChem.EmbedMolecule(mol_withh)
lines = Chem.MolToMolBlock(mol_withh)
lines = lines.split("\n")
headdata = lines[3].replace("\n", "").split(" ")
headdata = list(filter(None, headdata))
natoms = int(lines[3][0:3])
nbonds = int(lines[3][3:6])
if natoms > 99:
    for i in range(4 + natoms, 4 + natoms + nbonds):
        lines[i] = lines[i][:3] + ' ' + lines[i][3:]
file = open("cococ_1.sdf", "w")
file.write("\n".join(lines))
file.close()
