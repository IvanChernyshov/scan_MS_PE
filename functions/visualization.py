'''
Vizualization functions for MS data
'''

#%% Imports

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import colors, gridspec

from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage


#%% Functions

def interpolate_color(x):
    '''
    Interpolates x
    '''
    return np.interp(x = x, xp = [0, 255], fp = [0, 1])


def get_cdict(gradient):
    '''
    Transforms list of colors to cdict for matplotlib
    '''
    gradient = [[interpolate_color(_) for _ in color] for color in gradient]
    cdict = {}
    for i, channel in enumerate(['red', 'green', 'blue']):
        cdict[channel] = []
        for j, coord in enumerate(np.linspace(0, 1, len(gradient))):
            cdict[channel].append( [coord, gradient[j][i], gradient[j][i]] )
    
    return cdict


def plot_ms_heatmap(df, max_mz = 250, min_intensity = 100,
                    show_alkanes = False, **kwargs):
    '''
    Plots all mass spectra as heatmap
    '''
    # get cmap
    gradient = [[255, 255, 255], [169,  23,  69], [244, 109,  69],
                [253, 219, 127], [230, 241, 146], [112, 198, 162],
                [ 64,  57, 144]]
    cdict = get_cdict(gradient)
    cmap = colors.LinearSegmentedColormap('cmap', segmentdata = cdict)
    # prepare df
    X = df.sort_index(ascending = False)
    X.reset_index(inplace = True, drop = True)
    cols = np.array([m for m in df.columns[1:] if m <= max_mz])
    X = np.array(X.loc[:,cols])
    # get mass / time dependence
    ms = cols
    m_ints = np.sum(X, axis = 0)
    m_ints = [100*m/sum(m_ints) for m in m_ints]
    ts = df.index
    t_ints = np.sum(X, axis = 1)
    t_ints = [100*t/sum(t_ints) for t in t_ints]
    # final df preparations
    X[X < min_intensity] = min_intensity
    vmin, vmax = X.min(), X.max()
    if sum(df.temp.isna()):
        extent = (min(cols), max(cols), min(df.index), max(df.index))
    else:
        extent = (min(cols), max(cols), min(df.temp), max(df.temp))
    # initialize subplots
    fig = plt.figure(figsize = (6, 4.5), dpi = 120)
    gs = gridspec.GridSpec(2, 2, width_ratios = [4, 1], height_ratios = [1, 4])
    gs.update(wspace = 0.025, hspace = 0.05)
    # main plot
    ax = plt.subplot(gs[1, 0])
    main = ax.imshow(X, extent = extent,
                     norm = colors.LogNorm(vmin = vmin, vmax = vmax),
                     cmap = cmap, aspect = 'auto')
    # show alkanes
    if show_alkanes:
        picks = [14*n + 2 for n in range(2, int(max_mz/14) + 1)]
        for pick in picks:
            ax.axvline(x = pick, linestyle = '--', color = 'black', lw = 1)
    # mass / time distribution
    axt = plt.subplot(gs[1, 1], frameon = False, sharey = ax,
                      xticks = [])
    if sum(df.temp.isna()):
        axt.plot(t_ints, ts[::-1])
    else:
        axt.plot(t_ints, df.temp[::-1])
    plt.setp(axt.get_yticklabels(), visible = False)
    plt.setp(axt.get_yticklines(), visible = False)
    # fraction distribution
    axm = plt.subplot(gs[0, 0], frameon = False, sharex = ax,
                      yticks = [])
    axm.plot(ms, m_ints)
    plt.setp(axm.get_xticklabels(), visible = False)
    plt.setp(axm.get_xticklines(), visible = False)
    # colorbar
    plt.colorbar(main, ax = axt)
    
    return fig


def plot_fraction_types(ftypes, temperatures, step_size = 50, **kwargs):
    '''
    Plots fraction types
    '''
    # get percentage
    total = sum([sum(ftype) for ftype in ftypes.values()])
    ps = {key: sum(val)/total for key, val in ftypes.items()}
    # plot figure
    fig, ax = plt.subplots(figsize = (6, 4.5), dpi = 120)
    if temperatures is None:
        xs = list(range(0, step_size*len(list(ftypes.values())[0]), step_size))
    elif sum(temperatures.isna()):
        xs = list(range(0, step_size*len(list(ftypes.values())[0]), step_size))
    else:
        idxs = list(range(0, step_size*len(list(ftypes.values())[0]), step_size))
        xs = [t for i, t in enumerate(temperatures) if i in idxs]
    legend = []
    for ftype, intensities in ftypes.items():
        plt.plot(xs, intensities)
        legend.append(f'{ftype}: {100*ps[ftype]:.1f}%')
    plt.legend(legend, loc = 'best', fontsize = 'small')
    
    return fig


def plot_main_compounds(quantities, fracs, basis, temperatures,
                        step_size = 50):
    '''
    Plots dependence of compound release on time/temperature
    '''
    # set temperatures
    if temperatures is None:
        xs = list(range(0, step_size*fracs.shape[0], step_size))
    elif sum(temperatures.isna()):
        xs = list(range(0, step_size*fracs.shape[0], step_size))
    else:
        idxs = list(range(0, step_size*fracs.shape[0], step_size))
        xs = [t for i, t in enumerate(temperatures) if i in idxs]
    # set figures
    figs = {}
    for t in basis.type.drop_duplicates():
        fig, ax = plt.subplots(figsize = (6, 3.5), dpi = 120)
        names = quantities.name[quantities.name.isin(basis.name[basis.type == t])].iloc[:6]
        legend = []
        for name in names:
            idx = basis.index[basis.name == name].values[0]
            ax.plot(xs, fracs[:,idx])
            q = quantities.percentage[quantities.name == name].values[0]
            legend.append(f'{name}: {q:.2f}%')
        plt.legend(legend, loc = 'best', fontsize = 'x-small')
        figs[t] = fig
    
    return figs


def print_main_compounds(quantities, basis):
    '''
    Plots all compounds and their percentages; support function
    '''
    mols = [Chem.MolFromSmiles(basis.smiles[basis.name == name].values[0]) for name in quantities.name]
    legends = [f'{p:.2f} %' for p in quantities.percentage]
    fig = MolsToGridImage(mols, molsPerRow = 5, subImgSize = (180, 135), legends = legends, maxMols = 100)
    
    return fig


