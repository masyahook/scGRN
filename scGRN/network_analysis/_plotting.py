import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from ..utils import _GGPLOT_COLORS as colors


def plot_avail_cell_types(meta: pd.DataFrame, save_as: str = None):
    """
    Plot the cell type distribution across patients as a heatmap.

    :param meta: The metadata data as a Pandas dataframe
    :param save_as: The file path to saved file, None if do not save the plot
    """

    df = meta.T.iloc[3:, :]
    N_c, N_p = df.shape

    for col in df.columns:
        df[col] = df[col].astype('float')
    df = df.assign(sum_col=df.sum(1)).sort_values('sum_col', ascending=False).drop(columns=['sum_col'])

    with sns.axes_style("white"):
        g = sns.JointGrid(ratio=15)
        g.ax_marg_y.cla()
        g.ax_marg_x.cla()
        sns.heatmap(data=df, ax=g.ax_joint, cbar=False, vmin=0, annot=True,
                    fmt='.0f', cmap='crest', mask=df.isnull())

        bars_y = g.ax_marg_y.barh(np.arange(0.5, N_c), df.sum(axis=1), color='lightgray')
        bars_x = g.ax_marg_x.bar(np.arange(0.5, N_p), df.sum(axis=0), color='lightgray')

        for bar in bars_x:
            g.ax_marg_x.annotate(
                int(bar.get_height()),
                (bar.get_x() + bar.get_width() / 2, 0), ha='center', va='center',
                size=12, xytext=(0, 8),
                textcoords='offset points'
            )

        for bar in bars_y:
            g.ax_marg_y.annotate(
                int(bar.get_width()),
                (0, bar.get_y() + bar.get_height() / 2), ha='left', va='center',
                size=12, xytext=(4, 0),
                textcoords='offset points'
            )

        g.ax_joint.set_xticks(np.arange(0.5, N_p))
        g.ax_joint.set_xticklabels(df.columns, rotation=0)
        g.ax_joint.set_yticks(np.arange(0.5, N_c))
        g.ax_joint.set_yticklabels(df.index.map(lambda x: x.replace('_', ' ')), rotation=0)

        # remove ticks between heatmap and histograms
        g.ax_marg_x.tick_params(axis='x', bottom=False, labelbottom=False)
        g.ax_marg_y.tick_params(axis='y', left=False, labelleft=False)
        # remove ticks showing the heights of the histograms
        g.ax_marg_x.tick_params(axis='y', left=False, labelleft=False)
        g.ax_marg_y.tick_params(axis='x', bottom=False, labelbottom=False)

        g.fig.set_size_inches(16, 8)  # jointplot creates its own figure, the size can only be changed afterwards
        # g.fig.subplots_adjust(hspace=0.3) # optionally more space for the tick labels
        g.fig.subplots_adjust(hspace=0.1, wspace=0.05)  # less spaced needed when there are no tick labels

        for tick_label in g.ax_joint.get_xticklabels():
            if meta['group'][tick_label.get_text()] == 'S':
                tick_label.set_color(colors['red'])
            elif meta['group'][tick_label.get_text()] == 'M':
                tick_label.set_color(colors['yellow'])
            else:
                tick_label.set_color(colors['green'])

        g.ax_joint.set_xlabel('')

    if save_as is not None:
        plt.savefig(save_as, bbox_inches='tight')
