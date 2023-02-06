'''
Functions for qualitative composition analysis of MS
'''

#%% Imports

import pandas as pd
import numpy as np
from scipy.optimize import least_squares

from rdkit import Chem


#%% Functions

def read_ms_data(path, temps, times):
    '''
    Returns MS data as dataframe (columns are m/z)
    '''
    df = pd.read_csv(path, skiprows = 4)
    df.drop(columns = [df.columns[0], df.columns[-1]], inplace = True)
    df.columns = [int(_.strip()) for _ in df.columns]
    # add temperature scale
    df.insert(0, 'temp', np.nan)
    if (temps is not None) and (times is not None):
        idx_end = df.index[-1]
        idx_start = round(len(df.index)*times[0]/times[1])
        df = df.loc[idx_start:idx_end]
        df.reset_index(drop = True, inplace = True)
        df['temp'] = np.linspace(temps[0], temps[1], len(df))
    
    return df


def get_mol_type(smi):
    '''
    Classifies mol as alkane / alkene / polyene / aromatics
    '''
    mol = Chem.MolFromSmiles(smi)
    bs = [str(b.GetBondType()) for b in mol.GetBonds()]
    if len(bs) == 0:
        return 'alkanes'
    elif 'AROMATIC' in bs:
        return 'aromatics'
    elif 'DOUBLE' in bs:
        return 'alkenes' if bs.count('DOUBLE') == 1 else 'polyenes'
    elif 'SINGLE' in bs and len(set(bs)) == 1:
        return 'alkanes'
    
    return 'others'


def get_ms_basis(path_basis, norm_area = True, min_mz = 15):
    '''
    Prepares data for spectrum deconvolution
    '''
    basis = pd.read_csv(path_basis)
    # get max_mz
    max_mz = max([max([int(_.split(',')[0]) for _ in ms.split()]) for ms in basis.MS])
    specs = []
    for idx, ms in zip(basis.index, basis.MS):
        signals = [_.split(',') for _ in ms.split()]
        signals = {int(key): int(val) for key, val in signals if int(key) >= min_mz}
        spec = np.zeros(max_mz + 1)
        for mass, intens in signals.items():
            spec[mass] = intens
        spec = spec[min_mz:]
        if norm_area:
            spec *= 10000 / sum(spec)
        specs.append(spec)
    basis['MS'] = specs
    
    return basis


def deconvolute_ms(df, basis, step_size = 50):
    '''
    Decyphers MS data
    '''
    df = df.loc[range(0, len(df), step_size)]
    # prepare spec decyphering
    A = np.array(list(basis.MS)).transpose()
    n, m = A.shape
    x0 = np.random.rand(m)
    fracs = []
    for ii, idx in enumerate(df.index):
        b = np.array(df.loc[idx,df.columns[1:]])[:len(basis.MS[0])]
        x = least_squares(lambda x: A @ x - b, x0, bounds = (0, float('inf')))
        fracs.append(x.x)
    fracs = np.array(fracs)
    # get types
    types = {t: basis.index[basis.type == t] for t in basis.type.drop_duplicates()}
    ftypes = {ftype: fracs[:,idxs].sum(axis = 1) for ftype, idxs in types.items()}
    
    return fracs, ftypes


def get_main_compounds(fracs, basis, min_percentage = 0.1, print_data = False):
    '''
    Returns table containing main compounds
    '''
    quantity = []
    name = []
    for idx in range(fracs.shape[1]):
        quantity.append(sum(fracs[:,idx]))
        name.append(basis.name[idx])
    quantity = [round(100*q/sum(quantity), 2) for q in quantity]
    q = pd.DataFrame({'name': name, 'percentage': quantity})
    q = q.sort_values('percentage', ascending = False)
    q = q.loc[q.percentage >= min_percentage]
    q = q.reset_index(drop = True)
    if print_data:
        pd.options.display.max_rows = 160
        print(q)
    
    return q


