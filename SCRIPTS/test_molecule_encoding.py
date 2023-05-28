import unittest
from rdkit import Chem

from encode_molecule import *

class TestMoleculeEncoding(unittest.TestCase):

    def test_atom_encoders(self):
        self.assertEqual(encode_atoms_HCNO.n_features(), 5)
        self.assertEqual(encode_atoms.n_features(), 43)

    def test_sum(self):
        mol = Chem.MolFromSmiles("O=C=O")
        atom_features = mol_atoms_to_data(mol, encode_atoms_HCNO)
        self.assertEqual(atom_features[0].tolist(), [0.0, 0.0, 0.0, 1.0, 0.0]) # O
        self.assertEqual(atom_features[1].tolist(), [0.0, 1.0, 0.0, 0.0, 0.0]) # O
        self.assertEqual(atom_features[2].tolist(), [0.0, 0.0, 0.0, 1.0, 0.0]) # O
        self.assertEqual(len(atom_features), 3)


if __name__ == '__main__':
    unittest.main()
