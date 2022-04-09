import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from favs.plot import stripplot_with_median_ci
import dabest as db
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

order = ['fed', 'deprived']
matings_states = ['virgin', 'mated']
regions = {
                'motor': [19,44,73,74]
}

df = pd.read_csv('./peaks.csv').set_index('Unnamed: 0')
df['group'] = df['MS'] +' ' + df['FS'] ## grouping for DaBest

violinplot_kwargs = {'widths':0.5, 'vert':True, 'showextrema':False, 'showmedians':False}
title_colors = []
bootstrap_df = {'region': [],
                'state': [],
                'mean_diff': [],
                'ci_high': [],
                'ci_low': [],}
for each_panel in regions.keys():
    print(each_panel)
    #df = pd.read_csv('bootstrap_data_{}.csv'.format(each_panel))
    f,axes = plt.subplots(ncols=len(regions[each_panel]),figsize=(len(regions[each_panel]),3), sharey=True)
    for i, (ax, region) in enumerate(zip(axes,regions[each_panel])):
        #ax.set_frame_on(False)
        if i==0:
            ax.set_ylabel('Mean difference')
        ax.get_xaxis().set_visible(False)
        for tick in range(2):
            analysis_of_long_df = db.load(df.loc[(df.region==region)&(df.stimulus=='yeast')], idx=(("virgin fed", "mated fed"),
                                                                                                   ("virgin deprived", "mated deprived")),
                                                                                                      x="group", y="peak",resamples=10000)
            results = analysis_of_long_df.mean_diff.results
            current_bootstrap = results.bootstraps[tick]
            current_effsize   = results.difference[tick]
            current_ci_low    = results.bca_low[tick]
            current_ci_high   = results.bca_high[tick]
            bootstrap_df['region'].append(region)
            bootstrap_df['state'].append(order[tick])
            bootstrap_df['mean_diff'].append(current_effsize)
            bootstrap_df['ci_high'].append(current_ci_high)
            bootstrap_df['ci_low'].append(current_ci_low)
            v = ax.violinplot(current_bootstrap[~np.isinf(current_bootstrap)],
                                 positions=[tick],
                                 **violinplot_kwargs)
            if tick==0:
                fc='lightgray'
            else:
                fc='#9986A5'
            halfviolin(v, fill_color=fc, alpha=.8)

            # Plot the effect size.
            ax.plot([tick], current_effsize, marker='o',
                   color='k',
                   markersize=9)
            # Plot the confidence interval.
            ax.plot([tick, tick],
                   [current_ci_low, current_ci_high],
                   linestyle="-",
                   color='k',
                   linewidth=2)

            ax.axhline(y=0,xmin=0,xmax=1,ls='dashed',lw=.5, color='k',zorder=0)
            ax.set_ylim([-.22,0.62])
            ax.set_title(region)

    f.suptitle(each_panel,x=0.15,ha='left')
    f.tight_layout()
    f.savefig('figS13_boot_deprivation_{}.pdf'.format(each_panel), dpi=300)
    plt.close()
bootstrap_df = pd.DataFrame(bootstrap_df)
bootstrap_df.to_csv('figS13_bootstrap_results_deprivation.csv')
