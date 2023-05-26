# RDkit
from rdkit import Chem
from rdkit.Chem.rdmolops import GetAdjacencyMatrix

# load data as tsv with pandas
import pandas as pd

import numpy as np
import torch

from encoder import Encoder

encode_atoms = Encoder(*['C','N','O','S','F','Si','P','Cl','Br','Mg','Na','Ca',
		    'Fe','As','Al','I', 'B','V','K','Tl','Yb','Sb','Sn','Ag',
		    'Pd','Co','Se','Ti','Zn', 'Li','Ge','Cu','Au','Ni','Cd',
		    'In','Mn','Zr','Cr','Pt','Hg','Pb','Unknown'])

def mol_atoms_to_data(mol_object):

    n_atoms = mol_object.GetNumAtoms()
    n_atom_features = encode_atoms.n_features()
    X = np.zeros((n_atoms, n_atom_features))
    for atom in mol_object.GetAtoms():
        s = atom.GetSymbol()
        X[atom.GetIdx(), :] = np.array(encode_atoms.encode_binary(s))

    X = torch.tensor(X, dtype = torch.float)

    return X

def mol_bonds_to_data(mol_object):

# construct edge index array E of shape (2, n_edges)
    (rows, cols) = np.nonzero(GetAdjacencyMatrix(mol_object))


def smiles_to_graph_data(x_smiles):

    all_data = []
    for smiles in x_smiles:
        mol = Chem.MolFromSmiles(smiles)
        X = mol_atoms_to_data(mol)
        all_data.append(X)
    return all_data

#mol_atoms_to_data("CCN1C(=O)/C(=C2\SC(=S)N(CCCOC)C2=O)c2ccccc21")

df = pd.read_csv('../INPUTS/data/cyp2c19_veith.tab', sep='\t')
x_smiles = df['Drug']
y = df['Y']
x_features = smiles_to_graph_data(x_smiles[:1])

print(len(x_features))

