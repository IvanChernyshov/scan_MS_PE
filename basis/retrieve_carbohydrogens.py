'''
Getting MS spectra for CnHm systems from NIST Chemistry WebBook
'''

#%% Imports

import re
import pandas as pd
import nistchempy as nist


#%% Get compounds

# get bruttos
bruttos = []
for nC in range(1, 21):
    for nH in range(2, nC*2+3, 2):
        bruttos.append(f'C{nC}H{nH}')

# retrieving codes
search = nist.Search(Units = 'SI', MatchIso = True, NoIon = True, cMS = True)
IDs = set()
for brutto in bruttos:
    print(brutto)
    search.find_compounds(brutto, 'formula')
    if not search.success or search.lost:
        raise f'Bad brutto: {brutto}'
    IDs.update(search.IDs)


#%% Downloading mass spectra

# load info
info = []
for ii, ID in enumerate(IDs):
    if not ii % 10:
        print(ii)
    X = nist.Compound(ID)
    if 'cMS' not in X.data_refs:
        continue
    X.get_ms_spectra()
    if not X.MS:
        continue
    spec = X.MS[0]
    MS = re.search(r'(\d+,\d+\s+)+', spec.jdx_text).group(0).replace('\n', ' ').strip()
    info.append( [X.ID, X.name, X.inchi, X.mol_weight, MS] )

# get table
df = pd.DataFrame(info, columns = ['ID', 'name', 'inchi', 'MW', 'MS'])
df.to_csv('csvs/CnHm_raw.csv', index = None)


