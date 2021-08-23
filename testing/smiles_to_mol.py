from rdkit import Chem
from rdkit.Chem import AllChem

SMILES = "[H][C@]1(CCC2)[C@@]2([H])CCC1"

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
file = open("molfromsmiles.sdf", "w")
file.write("\n".join(lines))
file.close()
