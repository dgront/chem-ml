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

encode_atoms_HCNO = Encoder(*['H','C','N','O','Unknown'])

encode_bond_types = Encoder(*[Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE,
			    Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC])

def mol_atoms_to_data(mol_object, encoder=encode_atoms):

    n_atoms = mol_object.GetNumAtoms()
    n_atom_features = encoder.n_features()
    X = np.zeros((n_atoms, n_atom_features))
    for atom in mol_object.GetAtoms():
        s = atom.GetSymbol()
        X[atom.GetIdx(), :] = np.array(encoder.encode_binary(s))

    X = torch.tensor(X, dtype = torch.float)

    return X

def bond_features(a_bond):
    bond_features_vector = encode_bond_types.encode_binary(a_bond.GetBondType())
    return np.array(bond_features_vector)


def mol_bonds_to_data(mol_object):

    # --- construct edge index array E of shape (2, n_edges)
    (rows, cols) = np.nonzero(GetAdjacencyMatrix(mol_object))
    n_edges = len(rows)
    n_edge_features = encode_bond_types.n_features()
    edge_features = np.zeros((n_edges, n_edge_features))
    for (k, (i,j)) in enumerate(zip(rows, cols)):
        edge_features[k, :] = np.array(bond_features(mol_object.GetBondBetweenAtoms(int(i),int(j))))

    return (rows, cols, edge_features)

def molecule_to_graph_data(mol_object):

    X = mol_atoms_to_data(mol_object)
    EF = mol_bonds_to_data(mol_object)
    return (X, EF)

def smiles_to_graph_data(x_smilest):

    all_data = []
    for smiles in x_smiles:
        mol = Chem.MolFromSmiles(smiles)
        molecule_to_graph_data(mol)
    return all_data

if __name__ == "__main__":
    #mol_atoms_to_data("CCN1C(=O)/C(=C2\SC(=S)N(CCCOC)C2=O)c2ccccc21")

    df = pd.read_csv('../INPUTS/data/cyp2c19_veith.tab', sep='\t')
    x_smiles = df['Drug']
    y = df['Y']
    x_features = smiles_to_graph_data(x_smiles[:1])

    print(len(x_features))

