"""
LOADING DATA
"""
import pandas as pd
df = pd.read_csv('./data_onsets.csv')
df['state'] = df['MS'] + ' ' + df['FS'] ## grouping for DaBest
print(df.columns)
# region: index of region
# group: combination of mating and metabolic state (mated fed & mated deprived)
# fly: index of fly (N=19)
# stimulus: taste stimulus (only yeast for filtered onsets)
# concentration: tastant concentration (10% for yeast)
# stimulus_nr: index no of stimulus
# recording: recording name
# FS: Feeding State (fed & deprived)
# MS: Mating State (only mated for filtered onsets)
# onset: filtered onset values [s]
# cat: region categories (sensory, motor & other)
print(df.region.unique())

"""
DEFINITIONS
"""
order = ['fed', 'deprived']
Index(['region', 'state', 'fly', 'stimulus', 'concentration', 'stimulus_nr',
       'recording', 'FS', 'MS', 'onset', 'cat'],
      dtype='object')
[ 5  6 14 16 17 19 22 25 26 41 43 44 45 56 57 58 59 60 62 63 66 68 72 73
 74 77 78  4 24 30 37 39 53 67 70]
import dabest as db
from tqdm.notebook import tqdm

effect_sizes_defs = ['mean_diff', 'cohens_d']
dfs = {k:[] for k in effect_sizes_defs}
effect_sizes = {'region': [], 'cat': []}
for k in effect_sizes_defs:
    effect_sizes[k] = []
for region in tqdm(df.region.unique()):
    analysis_of_long_df = db.load(df.loc[df.region==region],
                                  idx=('mated fed', 'mated deprived'),
                                  x='state',
                                  y='onset',
                                  resamples=10000)
    effect_sizes['region'].append(region)
    effect_sizes['cat'].append(df.loc[df.region==region].cat.unique()[0])
    for key in effect_sizes_defs:
        val = float(getattr(analysis_of_long_df,key).results.difference.round(3))
        #print(region, key, val)
        effect_sizes[key].append(val)
        bsdata = getattr(analysis_of_long_df,key).results.bootstraps
        for i, row in bsdata.iteritems():
            _df = pd.DataFrame({'value': row})
            _df['key'] = '{}_{}'.format(region, i)
            dfs[key].append(_df)
for key in effect_sizes_defs:
    dfs[key] = pd.concat(dfs[key], ignore_index=True)
    dfs[key].to_csv('{}_bootstrap_data_onsets.csv'.format(key))
effect_sizes = pd.DataFrame(effect_sizes)
effect_sizes.to_csv('effect_sizes_onsets.csv'.format(key))

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from favs.plot import stripplot_with_median_ci
import sys

def halfviolin(v, half='right', fill_color='k', alpha=1,
                line_color='k', line_width=0):
    import numpy as np

    for b in v['bodies']:
        V = b.get_paths()[0].vertices

        mean_vertical = np.mean(V[:, 0])
        mean_horizontal = np.mean(V[:, 1])

        if half == 'right':
            V[:, 0] = np.clip(V[:, 0], mean_vertical, np.inf)
        elif half == 'left':
            V[:, 0] = np.clip(V[:, 0], -np.inf, mean_vertical)
        elif half == 'bottom':
            V[:, 1] = np.clip(V[:, 1], -np.inf, mean_horizontal)
        elif half == 'top':
            V[:, 1] = np.clip(V[:, 1], mean_horizontal, np.inf)

        b.set_color(fill_color)
        b.set_alpha(alpha)
        b.set_edgecolor(line_color)
        b.set_linewidth(line_width)
effect_size_type = 'mean_diff'
effect_sizes = pd.read_csv('effect_sizes_onsets.csv')
print(effect_sizes)
violinplot_kwargs = {'widths':0.1, 'vert':True, 'showextrema':False, 'showmedians':False}
title_colors = []
colors = {'sensory': '#FCC765', 'motor': '#49A7CC', 'other': '#8C8C8C'}
dfs = pd.read_csv('{}_bootstrap_data_onsets.csv'.format(effect_size_type))
order = effect_sizes.sort_values(by=effect_size_type).region.values[:2]
cats = effect_sizes.sort_values(by=effect_size_type).cat.values[:2]
f,axes = plt.subplots(ncols=len(order),figsize=(len(order),3), sharey=True)
for i, (ax, region,cat) in tqdm(enumerate(zip(axes,order,cats)), total=len(order)):
    #ax.set_frame_on(False)
    if i==0:
        ax.set_ylabel('Mean differences')#('Cohen\'s d')
    ax.get_xaxis().set_visible(False)
    analysis_of_long_df = db.load(df.loc[df.region==region],
                                  idx=("mated fed", "mated deprived"),
                                  x="state",
                                  y="onset",
                                  resamples=10000)
    results = analysis_of_long_df.mean_diff.results
    current_bootstrap = results.bootstraps[0]
    current_effsize   = results.difference[0]
    current_ci_low    = results.bca_low[0]
    current_ci_high   = results.bca_high[0]
    print(current_effsize, results.pvalue_welch[0], results.pvalue_mann_whitney[0])
    v = ax.violinplot(current_bootstrap[(~np.isinf(current_bootstrap))&(~np.isnan(current_bootstrap))],
                         positions=[0],
                         **violinplot_kwargs)
    for pc in v['bodies']:
        pc.set_facecolor(colors[cat])
        #halfviolin(v, fill_color=fc, alpha=.8)

    # Plot the effect size.
    ax.plot([0], current_effsize, marker='o',
           color='k',
           markersize=9)
    # Plot the confidence interval.
    ax.plot([0, 0],
           [current_ci_low, current_ci_high],
           linestyle="-",
           color='k',
           linewidth=2)

    ax.axhline(y=0,xmin=0,xmax=1,ls='dashed',lw=.5, color='k',zorder=0)
    ax.set_ylim([-6,6])
    ax.set_title(region)
    sns.despine(ax=ax, trim=True, bottom=True)

f.tight_layout()
f.savefig('{}.pdf'.format(effect_size_type), dpi=300)
plt.close()
