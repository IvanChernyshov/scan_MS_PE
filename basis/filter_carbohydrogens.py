'''
Getting MS spectra for CnHm systems from NIST Chemistry WebBook
'''

#%% Imports

import itertools
import pandas as pd
from rdkit import Chem


#%% Functions

def count_elems(mol):
    '''
    Returns tuple of number of carbons, hydrogens, and others
    '''
    nC = 0; nH = 0; nX = 0
    for a in mol.GetAtoms():
        sym = a.GetAtomicNum()
        if sym == 6:
            nC += 1
        elif sym == 1:
            nH += 1
        else:
            nX += 1
        nH += a.GetNumImplicitHs()
    
    return (nC, nH, nX)


def filter_aliphatics(mol):
    '''
    Returns True, if mol is non-aromatic, non-cyclic,
    have single/double bonds only, and low branching degree
    '''
    # check cycle
    if mol.GetRingInfo().NumRings():
        return False
    # check bonds
    bs = set([str(b.GetBondType()) for b in mol.GetBonds()])
    if bs.difference({'SINGLE', 'DOUBLE'}):
        return False
    # check branching
    nC = len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 6])
    max_dist = 0
    for i, j in itertools.combinations(range(mol.GetNumAtoms()), r = 2):
        dist = len(Chem.GetShortestPath(mol, i, j))
        if dist > max_dist:
            max_dist = dist
    if nC - max_dist > 2:
        return False
    
    return True


def filter_aromatics(mol):
    '''
    Returns True, if mol is aromatic, contains one aromatic block
    (biphenyl contains 2), and no more than 4 non-aromatic carbons
    '''
    if not mol.GetSubstructMatches(Chem.MolFromSmarts('c')):
        return False
    # check aliphatic carbons
    nCal = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C')))
    if nCal > 4:
        return False
    # check cycle blocks
    bs = [b.GetIdx() for b in mol.GetBonds() if not b.IsInRing()]
    if not bs:
        n_frags = int(bool(mol.GetRingInfo().NumRings()))
    else:
        frags = Chem.GetMolFrags(Chem.FragmentOnBonds(mol, bs), asMols = True)
        n_frags = len([f for f in frags if f.GetRingInfo().NumRings()])
    if n_frags > 1:
        return False
    
    return True


def get_mol_type(mol):
    '''
    Classifies mol as alkane / alkene / polyene / aromatics
    '''
    bs = [str(b.GetBondType()) for b in mol.GetBonds()]
    if len(bs) == 0: # methane
        return 'alkanes'
    elif 'AROMATIC' in bs:
        return 'aromatics'
    elif 'DOUBLE' in bs:
        return 'alkenes' if bs.count('DOUBLE') == 1 else 'polyenes'
    elif 'SINGLE' in bs and len(set(bs)) == 1:
        return 'alkanes'
    
    return 'others'



#%% Filter data

# read
df = pd.read_csv('csvs/CnHm_raw.csv')

# one mol for unique non-isomeric SMILES (with the shortest name)
df = df.loc[~df.inchi.isna()]
df['smiles'] = [Chem.MolToSmiles(Chem.MolFromInchi(_),
                                 isomericSmiles = False,
                                 canonical = True) for _ in df.inchi]
df['name_len'] = [len(_) for _ in df.name]
df = df.loc[df.groupby(['smiles']).name_len.idxmin()]

# drop bad brutto
df['nC'], df['nH'], df['nX'] = zip(*[count_elems(Chem.MolFromSmiles(_)) for _ in df.smiles])
df = df.loc[df.nX == 0]

# filter aromatics and aliphatics
mols = [Chem.MolFromSmiles(_) for _ in df.smiles]
idxs = [i for i, mol in enumerate(mols) if filter_aliphatics(mol) or filter_aromatics(mol)]
df = df.iloc[idxs]

# classify systems
mols = [Chem.MolFromSmiles(_) for _ in df.smiles]
df['type'] = [get_mol_type(_) for _ in mols]

# sort and save
df.sort_values(['type', 'nC', 'nH', 'name'], ascending = [True, True, False, True], inplace = True)
df.reset_index(drop = True, inplace = True)

# reorder columns and save
df = df.loc[:, ['name', 'type', 'smiles', 'nC', 'nH', 'MW', 'MS']]
df.to_csv('csvs/CnHm_filtered.csv')


