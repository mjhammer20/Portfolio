# Imports
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw

# Class Def

class rdkit_essentials():
    '''Class housing essential RDKit functions for basic cheminformatics tasks.'''

    def __init__(self):
        self.ms = []
        self.fps = []

    def mol_from_smiles(self, smiles):
        '''Generate RDKit Mol objects from SMILES strings.'''
        self.ms = [Chem.MolFromSmiles(smile) for smile in smiles]
    
    def generate_fingerprints(self):
        '''Generate fingerprints for the molecules.'''
        fpgen = AllChem.GetRDKitFPGenerator()
        self.fps = [fpgen.GetFingerprint(m) for m in self.ms]
    
    def calculate_similarity(self, index1, index2):
        '''Calculate Tanimoto similarity between two fingerprints.'''
        if index1 < len(self.fps) and index2 < len(self.fps):
            return DataStructs.TanimotoSimilarity(self.fps[index1], self.fps[index2])
        else:
            raise IndexError("Index out of range for fingerprints.")
        
    def substructure_search(self, index1, index2):
        '''Check if one molecule contains a substructure of another.'''
        if index1 < len(self.ms) and index2 < len(self.ms):
            return self.ms[index1].HasSubstructMatch(self.ms[index2])
        else:
            raise IndexError("Index out of range for molecules.")
        
    def visualize_molecules(self, smiles):
        '''Visualize a molecule in grid.'''
        return Draw.MolsToGridImage(self.ms, molsPerRow=3, subImgSize=(200, 200), legends=[smiles[i] for i in range(len(self.ms))])
        
