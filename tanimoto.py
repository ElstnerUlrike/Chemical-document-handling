from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
import re

def calculate_tanimoto_similarity(smiles1, smiles2):
    # Convert SMILES to RDKit molecules
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        return None  # Return None if SMILES strings are invalid
    
    # Generate Morgan fingerprints (ECFP4)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, 2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, 2048)
    
    # Calculate Tanimoto similarity
    tanimoto_sim = DataStructs.TanimotoSimilarity(fp1, fp2)
    
    return tanimoto_sim
