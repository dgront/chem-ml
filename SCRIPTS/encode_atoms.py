# RDkit
from rdkit import Chem
from rdkit.Chem.rdmolops import GetAdjacencyMatrix

from encoder import Encoder

encode_atoms = Encoder(*['C','N','O','S','F','Si','P','Cl','Br','Mg','Na','Ca',
		    'Fe','As','Al','I', 'B','V','K','Tl','Yb','Sb','Sn','Ag',
		    'Pd','Co','Se','Ti','Zn', 'Li','Ge','Cu','Au','Ni','Cd',
		    'In','Mn','Zr','Cr','Pt','Hg','Pb','Unknown'])

def mol_atoms_to_data(x_smiles, y):

    data_list = []
    
    for (smiles, y_val) in zip(x_smiles, y):
        # convert SMILES to RDKit mol object
        mol = Chem.MolFromSmiles(smiles)

        # get feature dimensions
        n_atoms = mol.GetNumAtoms()
        atom_names = [atom.GetSymbol() for atom in mol.GetAtoms()]
        print(atom_names)
        enc = encode_atoms.encode_binary(*atom_names)
        enc = encode_atoms.encode_values(*atom_names)
        print(enc)

mol_atoms_to_data(["CCN1C(=O)/C(=C2\SC(=S)N(CCCOC)C2=O)c2ccccc21"],[1])