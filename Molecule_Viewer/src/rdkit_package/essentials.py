# Imports
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw

# Class Def

class rdkit_essentials():
    '''
    Class housing essential RDKit functions for basic cheminformatics tasks

    '''

    def __init__(self):
        self.smiles = []
        self.ms = []
        self.fps = []

    def import_smiles(self, smiles: str) -> list:
        '''
        Check SMILES for valid format. If valid, add to list and convert to RDKit Mol object.
        If invalid, raise ValueError.
        
        :input smiles: SMILES string.
        :type smiles: list[str]

        '''
        ms = Chem.MolFromSmiles(smiles)

        if not ms:
            raise ValueError("Invalid SMILES string provided.")
        self.smiles.append(smiles)
        self.ms.append(ms)

        return self.smiles
    
    def generate_fingerprints(self):
        '''
        Generate fingerprints for the molecules.
        
        '''

        fpgen = AllChem.GetRDKitFPGenerator()
        self.fps = [fpgen.GetFingerprint(m) for m in self.ms]
    
    def calculate_similarity(self, index1: int, index2: int):
        '''
        Calculate Tanimoto similarity between two fingerprints.
        
        :param index1: Index of the first molecule.
        :type index1: int
        :param index2: Index of the second molecule.
        :type index2: int
        :return: Tanimoto similarity score between two fingerprints.
        :rtype: float

        '''

        if index1 < len(self.fps) and index2 < len(self.fps):
            return DataStructs.TanimotoSimilarity(self.fps[index1], self.fps[index2])
        else:
            raise IndexError("Index out of range for fingerprints.")
        
    def substructure_search(self, index, substructure):
        '''
        Check if one molecule contains a substructure of another.
        
        :param index: Index of the molecule to search for substructure.
        :type index1: int
        :param substructure: SMILES string for substructure to search for.
        :type substructure: str
        :return: True if the first molecule contains a substructure of the second, False otherwise
        :rtype: bool

        '''
    
        sub_mol = Chem.MolFromSmiles(substructure)
        if not sub_mol:
            raise ValueError("Invalid substructure SMILES string provided.")
        has_substructure = self.ms[index].HasSubstructMatch(sub_mol)
        if not has_substructure:
            raise ValueError("Substructure not found in the molecule.")
        
        hit_atts = list(self.ms[index].GetSubstructMatch(sub_mol))
        hit_bonds = []
        for bond in sub_mol.GetBonds():
            bond_obj = self.ms[index].GetBondBetweenAtoms(hit_atts[bond.GetBeginAtomIdx()], hit_atts[bond.GetEndAtomIdx()])
            hit_bonds.append(bond_obj.GetIdx())
        img = Draw.rdMolDraw2D.MolDraw2DSVG(400, 200)
        
        return Draw.rdMolDraw2D.PrepareAndDrawMolecule(img, self.ms[index], highlightAtoms=hit_atts, highlightBonds=hit_bonds)
        
        
    def visualize_molecules(self):
        '''
        Generate image to visualize molecules in a grid.
        
        :param smiles: List of SMILES strings for the molecules.
        :type smiles: list[str]
        return: Image of the molecules in a grid format.
        :rtype: PIL.Image or None

        '''
        return Draw.MolsToGridImage(self.ms, molsPerRow=3, subImgSize=(200, 200), legends=[self.smiles[i] for i in range(len(self.ms))])
        
