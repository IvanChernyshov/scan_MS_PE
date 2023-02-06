'''
Getting MS spectra for CnHm systems from NIST Chemistry WebBook
'''

#%% Imports

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist#, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt


#%% Functions

def signals_to_spectrum(ms, max_mz):
    '''
    Returns np array of ms
    '''
    signals = [_.split(',') for _ in ms.split()]
    signals = {int(key): int(val) for key, val in signals}
    spec = np.zeros(max_mz + 1)
    for mass, intens in signals.items():
        spec[mass] = intens
    # normalize to one integral value
    spec *= 10000 / sum(spec)
    
    return spec



#%% Plot dendrograms

# read
df = pd.read_csv('csvs/CnHm_filtered.csv')

# clustering
for subtype in set(df.type):
    sub = df.loc[df.type == subtype]
    max_mz = int(100*(max(sub.MW) // 100 + 1))
    MSs = np.array([signals_to_spectrum(ms, max_mz) for ms in sub.MS])
    X = pdist(MSs)
    #plt.imshow(squareform(X))
    # make dendrogram
    plt.figure(figsize=(6, 24*len(sub)/271), dpi = 300)
    linked = linkage(X, 'single')
    dendrogram(linked, orientation = 'right', labels = list(sub.name),
               distance_sort = 'descending', show_leaf_counts = True)
    plt.subplots_adjust(left = 0.3)
    plt.savefig(f'dendrograms/{subtype}.png')


