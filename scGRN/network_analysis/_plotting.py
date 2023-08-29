"""Plotting functions for network analysis."""

from itertools import chain  # for aggregate functions
from typing import Dict, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap
from matplotlib.transforms import Bbox
from tqdm.notebook import tqdm
from wordcloud import WordCloud

from ..config import _ALPHA
from ..config import _GGPLOT_COLORS as colors
from ..config import _NODE_SIZE, _STOPWORDS
from ..utils import get_elipsis_mask, scale, scale_int
from ._auxiliary_data import load_gene_func_db


def plot_avail_cell_types(meta: pd.DataFrame, save_as: str = None):
    """
    Plot the cell type count distribution across patients as a heatmap.

    :param meta: The metadata data as a Pandas dataframe
    :param save_as: The file path to saved file, None if do not save the plot
    """

    df = meta.T.iloc[3:, :]
    N_c, N_p = df.shape

    for col in df.columns:
        df[col] = df[col].astype("float")
    df = (
        df.assign(sum_col=df.sum(1))
        .sort_values("sum_col", ascending=False)
        .drop(columns=["sum_col"])
    )

    with sns.axes_style("white"):
        g = sns.JointGrid(ratio=15)
        g.ax_marg_y.cla()
        g.ax_marg_x.cla()
        sns.heatmap(
            data=df,
            ax=g.ax_joint,
            cbar=False,
            vmin=0,
            annot=True,
            fmt=".0f",
            cmap="crest",
            mask=df.isnull(),
        )

        bars_y = g.ax_marg_y.barh(
            np.arange(0.5, N_c), df.sum(axis=1), color="lightgray"
        )
        bars_x = g.ax_marg_x.bar(np.arange(0.5, N_p), df.sum(axis=0), color="lightgray")

        for bar in bars_x:
            g.ax_marg_x.annotate(
                int(bar.get_height()),
                (bar.get_x() + bar.get_width() / 2, 0),
                ha="center",
                va="center",
                size=12,
                xytext=(0, 8),
                textcoords="offset points",
            )

        for bar in bars_y:
            g.ax_marg_y.annotate(
                int(bar.get_width()),
                (0, bar.get_y() + bar.get_height() / 2),
                ha="left",
                va="center",
                size=12,
                xytext=(4, 0),
                textcoords="offset points",
            )

        g.ax_joint.set_xticks(np.arange(0.5, N_p))
        g.ax_joint.set_xticklabels(df.columns, rotation=0)
        g.ax_joint.set_yticks(np.arange(0.5, N_c))
        g.ax_joint.set_yticklabels(
            df.index.map(lambda x: x.replace("_", " ")), rotation=0
        )

        # remove ticks between heatmap and histograms
        g.ax_marg_x.tick_params(axis="x", bottom=False, labelbottom=False)
        g.ax_marg_y.tick_params(axis="y", left=False, labelleft=False)
        # remove ticks showing the heights of the histograms
        g.ax_marg_x.tick_params(axis="y", left=False, labelleft=False)
        g.ax_marg_y.tick_params(axis="x", bottom=False, labelbottom=False)

        g.fig.set_size_inches(
            16, 8
        )  # jointplot creates its own figure, the size can only be changed afterwards
        # g.fig.subplots_adjust(hspace=0.3) # optionally more space for the tick labels
        g.fig.subplots_adjust(
            hspace=0.1, wspace=0.05
        )  # less spaced needed when there are no tick labels

        for tick_label in g.ax_joint.get_xticklabels():
            if meta["group"][tick_label.get_text()] == "S":
                tick_label.set_color(colors["red"])
            elif meta["group"][tick_label.get_text()] == "M":
                tick_label.set_color(colors["yellow"])
            else:
                tick_label.set_color(colors["green"])

        g.ax_joint.set_xlabel("")

    if save_as is not None:
        plt.savefig(save_as, bbox_inches="tight")


def fancy_draw_network_edge_labels(
    G: nx.DiGraph,
    pos: dict,
    edge_labels: dict = None,
    label_pos: float = 0.5,
    font_size: int = 10,
    font_color: str = "k",
    font_family: str = "sans-serif",
    font_weight: str = "normal",
    alpha: float = None,
    bbox: Bbox = None,
    horizontalalignment: str = "center",
    verticalalignment: str = "center",
    ax: Axes = None,
    rotate: bool = True,
    clip_on: bool = True,
    rad: float = 0,
) -> dict:
    """
    Draw edge labels with additional feature to make them curved.

    :param G: A NetworkX graph
    :param pos: A dictionary with nodes as keys and positions as values
    :param edge_labels: Edge labels in a dictionary of labels keyed by edge two-tuple
    :param label_pos: Position of edge label along edge (0=head, 0.5=center, 1=tail)
    :param font_size: Font size for text labels
    :param font_color: Font color string
    :param font_weight: Font weight
    :param font_family: Font family
    :param alpha: The text transparency
    :param bbox: Specify text box properties (e.g. shape, color etc.) for edge labels
    :param horizontalalignment: Horizontal alignment {'center', 'right', 'left'}
    :param verticalalignment: Vertical alignment {'center', 'top', 'bottom', 'baseline', 
        'center_baseline'}
    :param ax: Draw the graph in the specified Matplotlib axes
    :param rotate: Rotate edge labels to lie parallel to edges
    :param clip_on: Turn on clipping of edge labels at axis boundaries

    :return Dictionary of labels keyed by edge
    """

    if ax is None:
        ax = plt.gca()
    if edge_labels is None:
        labels = {(u, v): d for u, v, d in G.edges(data=True)}
    else:
        labels = edge_labels
    text_items = {}
    for (n1, n2), label in labels.items():
        (x1, y1) = pos[n1]
        (x2, y2) = pos[n2]
        (x, y) = (
            x1 * label_pos + x2 * (1.0 - label_pos),
            y1 * label_pos + y2 * (1.0 - label_pos),
        )
        pos_1 = ax.transData.transform(np.array(pos[n1]))
        pos_2 = ax.transData.transform(np.array(pos[n2]))
        linear_mid = 0.5 * pos_1 + 0.5 * pos_2
        d_pos = pos_2 - pos_1
        rotation_matrix = np.array([(0, 1), (-1, 0)])
        ctrl_1 = linear_mid + rad * rotation_matrix @ d_pos
        ctrl_mid_1 = 0.5 * pos_1 + 0.5 * ctrl_1
        ctrl_mid_2 = 0.5 * pos_2 + 0.5 * ctrl_1
        bezier_mid = 0.5 * ctrl_mid_1 + 0.5 * ctrl_mid_2
        (x, y) = ax.transData.inverted().transform(bezier_mid)

        if rotate:
            # in degrees
            angle = np.arctan2(y2 - y1, x2 - x1) / (2.0 * np.pi) * 360
            # make label orientation "right-side-up"
            if angle > 90:
                angle -= 180
            if angle < -90:
                angle += 180
            # transform data coordinate angle to screen coordinate angle
            xy = np.array((x, y))
            trans_angle = ax.transData.transform_angles(
                np.array((angle,)), xy.reshape((1, 2))
            )[0]
        else:
            trans_angle = 0.0
        # use default box of white with white border
        if bbox is None:
            bbox = dict(boxstyle="round", ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0))
        if not isinstance(label, str):
            label = str(label)  # this makes "1" and 1 labeled the same

        t = ax.text(
            x,
            y,
            label,
            size=font_size,
            color=font_color,
            family=font_family,
            weight=font_weight,
            alpha=alpha,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            rotation=trans_angle,
            transform=ax.transData,
            bbox=bbox,
            zorder=1,
            clip_on=clip_on,
        )
        text_items[(n1, n2)] = t

    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False,
    )

    return text_items


def draw_graph(
    G: nx.DiGraph,
    pos: dict,
    ax: Axes,
    TF_names: list = None,
    label_edges: bool = True,
    node_size: float = _NODE_SIZE,
    alpha: float = _ALPHA,
    if_alpha_edges: bool = False,
    plot_cmap: bool = True,
    cmap: ListedColormap = plt.cm.plasma,
    label_font_size: float = 12,
) -> Axes:
    """
    Draw gene regulatory network using NetworkX.

    :param G: A NetworkX graph
    :param pos: A dictionary with nodes as keys and positions as values
    :param ax: Draw the graph in the specified Matplotlib axes
    :param TF_names: A list of transcription factors to highlight on the plot
    :param label_edges: True if label edges with importance feature, False otherwise
    :param node_size: The size of the gene nodes
    :param alpha: The level of transparency (between 0 and 1)
    :param if_alpha_edges: True if make edges transparent, False otherwise
    :param plot_cmap: True if plot colormap indicating Spearman correlation
    :param cmap: The colormap to use
    :param label_font_size: The text font size

    :return The matplotlib axes with drawn graph
    """

    # Define colors
    blue = np.expand_dims(np.array([221, 232, 250]) / 256, axis=0)
    dark_blue = np.expand_dims(np.array([115, 141, 187]) / 256, axis=0)
    yellow = np.expand_dims(np.array([253, 242, 208]) / 256, axis=0)
    dark_yellow = np.expand_dims(np.array([209, 183, 101]) / 256, axis=0)

    # Drawing nodes
    nx.draw_networkx_nodes(
        G, pos, node_color=yellow, edgecolors=dark_yellow, ax=ax, node_size=node_size
    )
    if TF_names is not None:
        nx.draw_networkx_nodes(
            G.subgraph(TF_names),
            pos=pos,
            node_color=blue,
            edgecolors=dark_blue,
            ax=ax,
            node_size=node_size,
        )

    # Drawing node labels
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=label_font_size)

    # Drawing edges
    if G.edges():
        edges, importances = zip(*nx.get_edge_attributes(G, "importance").items())
        edges, rhos = zip(*nx.get_edge_attributes(G, "rho").items())
        widths = scale(importances, 1, 10)

        curved_mask = [reversed(edge) in G.edges() for edge in G.edges()]

        edge_meta = {
            "curved_edge": [e for m, e in zip(curved_mask, edges) if m],
            "curved_importance": [e for m, e in zip(curved_mask, importances) if m],
            "curved_rho": [e for m, e in zip(curved_mask, rhos) if m],
            "curved_width": [e for m, e in zip(curved_mask, widths) if m],
            "straight_edge": [e for m, e in zip(curved_mask, edges) if not m],
            "straight_importance": [
                e for m, e in zip(curved_mask, importances) if not m
            ],
            "straight_rho": [e for m, e in zip(curved_mask, rhos) if not m],
            "straight_width": [e for m, e in zip(curved_mask, widths) if not m],
        }

        mpl_straight_edges = nx.draw_networkx_edges(
            G,
            pos,
            ax=ax,
            edgelist=edge_meta["straight_edge"],
            arrowstyle="->",
            arrowsize=30,
            edge_color=edge_meta["straight_rho"],
            edge_cmap=cmap,
            width=edge_meta["straight_width"],
            node_size=node_size,
        )
        mpl_curved_edges = nx.draw_networkx_edges(
            G,
            pos,
            ax=ax,
            edgelist=edge_meta["curved_edge"],
            connectionstyle=f"arc3, rad = 0.25",
            arrowstyle="->",
            arrowsize=30,
            edge_color=edge_meta["curved_rho"],
            edge_cmap=cmap,
            width=edge_meta["curved_width"],
            node_size=node_size,
        )

        if mpl_curved_edges is None:
            mpl_curved_edges = []
        if mpl_straight_edges is None:
            mpl_straight_edges = []

        if plot_cmap:
            pc = mpl.collections.PatchCollection(
                mpl_straight_edges + mpl_curved_edges, cmap=cmap
            )
            pc.set_array(rhos)
            cbar = plt.colorbar(pc, shrink=0.5)
            cbar.ax.set_ylabel(
                "Spearman correlation", rotation=270, fontsize=17, labelpad=17
            )

        if label_edges:
            edge_weights = nx.get_edge_attributes(G, "importance")
            curved_edge_labels = {
                edge: f"{edge_weights[edge]:.1f}" for edge in edge_meta["curved_edge"]
            }
            straight_edge_labels = {
                edge: f"{edge_weights[edge]:.1f}" for edge in edge_meta["straight_edge"]
            }
            fancy_draw_network_edge_labels(
                G, pos, ax=ax, edge_labels=curved_edge_labels, rotate=False, rad=0.25
            )
            fancy_draw_network_edge_labels(
                G, pos, ax=ax, edge_labels=straight_edge_labels, rotate=False
            )

        if if_alpha_edges:
            edge_importances = [
                w_dict["importance"] for st, end, w_dict in G.edges(data=True)
            ]
            alpha_importances = [
                (w - min(edge_importances))
                / (max(edge_importances) - min(edge_importances))
                * alpha
                + alpha
                for w in edge_importances
            ]
            # set alpha value for each edge
            for i in range(len(alpha_importances)):
                edges[i].set_alpha(alpha_importances[i])

    plt.axis("off")

    return ax


def graph_stats_vs_num_cells(graph_stats: dict, save_as: str = None) -> np.ndarray:
    """
    Plot graph properties of inferred GRNs against the number of cells of corresponding scRNA-seq 
    matrix.

    :param graph_stats: A dictionary containing graph properties, look in _data_processing.py in 
        get_graph_stats() for details
    :param save_as: Name of the figure file, or None if not to save

    :return The numpy array of matplotlib axis
    """

    fontsize = 25

    with sns.axes_style("white"):
        f, ax = plt.subplots(1, 3, figsize=(33, 7))

        # Plotting num_nodes vs num_edges vs num_cells
        ax[0].scatter(
            graph_stats["all"]["num_cells"],
            graph_stats["all"]["num_nodes"],
            color=colors["green"],
            label="Gene-gene networks",
        )
        ax[0].scatter(
            graph_stats["ctx"]["num_cells"],
            graph_stats["ctx"]["num_nodes"],
            color=colors["red"],
            label="TF regulon networks",
        )

        ax1 = ax[0].twinx()
        ax1.scatter(
            graph_stats["all"]["num_cells"],
            graph_stats["all"]["num_edges"],
            color=colors["green"],
        )
        ax1.scatter(
            graph_stats["ctx"]["num_cells"],
            graph_stats["ctx"]["num_edges"],
            color=colors["red"],
        )

        ax[0].set_xlabel("Number of cells in data", fontsize=fontsize)
        ax[0].set_ylabel("Number of nodes in network", fontsize=fontsize)
        ax1.set_ylabel(
            "Number of edges in network", rotation=270, labelpad=20, fontsize=fontsize
        )

        ax[0].set_yscale("log")
        ax1.set_yscale("log")

        ax[0].grid(True)
        ax1.grid(True)

        # Plotting radius vs diameter vs num_cells
        ax[1].scatter(
            graph_stats["all"]["num_cells"],
            graph_stats["all"]["diameter"],
            color=colors["green"],
        )
        ax[1].scatter(
            graph_stats["ctx"]["num_cells"],
            graph_stats["ctx"]["diameter"],
            color=colors["red"],
        )

        ax[1].set_xlabel("Number of cells in data", fontsize=fontsize)
        ax[1].set_ylabel("Network diameter", fontsize=fontsize)

        ax[1].grid(True)

        # Plotting average_degree vs average_path_length vs num_cells
        ax[2].scatter(
            graph_stats["all"]["num_cells"],
            graph_stats["all"]["average_degree"],
            color=colors["green"],
        )
        ax[2].scatter(
            graph_stats["ctx"]["num_cells"],
            graph_stats["ctx"]["average_degree"],
            color=colors["red"],
        )

        ax1 = ax[2].twinx()
        ax1.scatter(
            graph_stats["all"]["num_cells"],
            graph_stats["all"]["average_path_length"],
            color=colors["green"],
        )
        ax1.scatter(
            graph_stats["ctx"]["num_cells"],
            graph_stats["ctx"]["average_path_length"],
            color=colors["red"],
        )

        ax[2].set_xlabel("Number of cells in data", fontsize=fontsize)
        ax[2].set_ylabel("Node average degree in network", fontsize=fontsize)
        ax1.set_ylabel(
            "Average path length in network",
            rotation=270,
            labelpad=25,
            fontsize=fontsize,
        )

        ax[2].grid(True)
        ax1.grid(True)

        handles, labels = ax[0].get_legend_handles_labels()
        f.legend(
            handles,
            labels,
            loc="upper center",
            ncol=2,
            prop={"size": fontsize},
            bbox_to_anchor=(0.5, 1.1),
        )

        for i in range(3):
            ax[i].set_xscale("log")

        plt.tight_layout()

        if save_as is not None:
            plt.savefig(save_as, bbox_inches="tight")

    return ax


def graph_edge_stats_vs_num_cells(graph_stats: dict, save_as: str = None) -> np.ndarray:
    """
    Plot graph edges properties of inferred GRNs against the number of cells of corresponding 
    scRNA-seq matrix.

    :param graph_stats: A dictionary containing graph properties, look in _data_processing.py in 
        get_graph_stats() for details
    :param save_as: Name of the figure file, or None if not to save

    :return The numpy array of matplotlib axis
    """

    fontsize = 20

    with sns.axes_style("white"):
        f, ax = plt.subplots(1, 3, figsize=(30, 7))

        # Plotting importance and rho
        importance_all, rho_all = np.random.choice(
            graph_stats["all"]["importances"], size=10000
        ), np.random.choice(graph_stats["all"]["rhos"], size=10000)
        importance_ctx, rho_ctx = np.random.choice(
            graph_stats["ctx"]["importances"], size=10000
        ), np.random.choice(graph_stats["ctx"]["rhos"], size=10000)
        ax[0].scatter(importance_all, rho_all, color=colors["green"], alpha=0.1)
        ax[0].scatter(importance_ctx, rho_ctx, color=colors["red"], alpha=0.1)

        ax[0].set_xlabel("Link importance", fontsize=fontsize)
        ax[0].set_ylabel("Link Spearman correlation", fontsize=fontsize)

        ax[0].grid(True)

        # Plotting importance median and std
        ax[1].scatter(
            graph_stats["all"]["num_cells"],
            graph_stats["all"]["median_importance"],
            color=colors["green"],
            label="Gene-gene networks",
        )
        ax[1].scatter(
            graph_stats["ctx"]["num_cells"],
            graph_stats["ctx"]["median_importance"],
            color=colors["red"],
            label="TF regulon networks",
        )

        ax[1].set_xlabel("Number of cells in data", fontsize=fontsize)
        ax[1].set_ylabel("Link importance median", fontsize=fontsize)

        ax[1].grid(True)

        # Plotting rho median and std
        ax[2].scatter(
            graph_stats["all"]["num_cells"],
            graph_stats["all"]["median_rho"],
            color=colors["green"],
        )
        ax[2].scatter(
            graph_stats["ctx"]["num_cells"],
            graph_stats["ctx"]["median_rho"],
            color=colors["red"],
        )

        ax[2].set_xlabel("Number of cells in data", fontsize=fontsize)
        ax[2].set_ylabel("Link Spearman correlation median", fontsize=fontsize)

        ax[2].grid(True)

        handles, labels = ax[1].get_legend_handles_labels()
        f.legend(handles, labels, loc="upper center", ncol=2, prop={"size": fontsize})

        if save_as is not None:
            plt.savefig(save_as, bbox_inches="tight")

    return ax


def graph_num_regulons_vs_num_cells(
    graph_stats: dict, save_as: str = None
) -> sns.JointGrid:
    """
    Plot the distribution of the number of inferred regulons against the number of cells of 
    corresponding scRNA-seq matrix.

    :param graph_stats: A dictionary containing graph properties, look in _data_processing.py in 
        get_graph_stats() for details
    :param save_as: Name of the figure file, or None if not to save

    :return The Seaborn JointGrid containing figure
    """

    with sns.axes_style("white"):
        # Plotting importance median and std
        g = sns.jointplot(
            x=graph_stats["ctx"]["num_cells"],
            y=graph_stats["ctx"]["num_tfs"],
            color=colors["blue"],
        )

        g.ax_marg_x.cla()
        g.ax_marg_x.axis("off")
        g.fig.set_size_inches(16, 8)

        dist_median = np.median(graph_stats["ctx"]["num_tfs"])
        g.ax_marg_y.axhline(dist_median, color="k", linestyle="dashed", linewidth=1)
        min_xlim, max_xlim = g.ax_marg_y.get_xlim()
        g.ax_marg_y.text(
            max_xlim * 0.95,
            dist_median * 0.95,
            "Median: {:.1f}".format(dist_median),
            horizontalalignment="right",
            verticalalignment="center",
            fontsize=15,
        )

        g.ax_joint.set_xlabel("Number of cells in data", fontsize=20)
        g.ax_joint.set_ylabel("Number of regulons", fontsize=20)

        g.ax_joint.grid(True)

        g.fig.suptitle("The number of inferred regulons", fontsize=22)

        if save_as is not None:
            plt.savefig(save_as, bbox_inches="tight")

    return g


def plot_cloud(
    G: nx.DiGraph,
    partition: Dict[str, int],
    squeezed_pos: Dict[str, Tuple[float, float]],
    ax: Axes,
    anno_db: str,
    display_func: bool = False,
    central_genes: Union[bool, Dict[int, Dict[str, float]]] = True,
    if_betweenness: bool = True,
    limit_anno_until: int = 50,
    k: int = 3000,
) -> Axes:
    """
    Plot word clouds depicting communities in the graph. Communities will be laid out on the
    periphery, with each node being a gene. Either gene list or gene function will be displayed on
    top of each community (i.e. cloud of nodes). If gene list is displayed `display_func=False`,
    then the size of the gene corresponds to the centrality of the gene inside corresponding
    community. If gene function is displayed `display_func=True`, then the size of the function
    corresponds to the frequency of the functional term attributed to the top central genes.

    :param G: NetworkX graph
    :param partition: A dictionary where key is the node name and value is the community number,
        same as `node_to_community`
    :param squeezed_pos: A dictionary containing node coordinates (preferably squeezed version to
        speed up computation, look for `squeeze_graph` function)
    :param ax: an matplotlib axis object
    :param anno_db: a database tag, could be MSigDB, GO, DoRothEA, or etc. Look in
        `load_gene_func_db` for the full list
    :param display_func: True if display the gene functions in the word cloud, False if display gene
        names
    :param central_genes: Could be a boolean or a dictionary of dictionaries.
        If boolean:
            True - use only top `limit_anno_until` important genes (based on centrality) for
                plotting
            False - otherwise
        If dictionary:
            Pass the central genes and scores for each community for plotting. The dictionary has
                the following format:
                central_genes = {
                    0: {
                        gene_name_1: gene_centrality_score_1,
                        gene_name_2: gene_centrality_score_2,
                        ...
                    },
                    ...
                }
    :param if_betweenness: True if use betweenness centrality as node importance score, False if use
        closeness centrality
    :param limit_anno_until: Number of genes to use to calculate wordcloud
    :param k: Use k nodes to estimate centrality

    :returns: The axis with the plotted communities as word clouds
    """

    # Loading the gene functional annotation
    gene_func = load_gene_func_db(anno_db, reload=False, as_series=True)

    # Reverting partition dict -> {group_1: [gene_1, gene_2, ...], group_2: [gene_3, ...], ...}
    partition_genes_ = {}
    for gene, i in partition.items():
        if i not in partition_genes_.keys():
            partition_genes_[i] = [gene]
        else:
            partition_genes_[i] += [gene]

    # Whether to filter the genes on which we compute the word cloud (most important genes)
    if central_genes is True:
        compute_centrality = (
            nx.betweenness_centrality if if_betweenness else nx.closeness_centrality
        )
        kwargs = {"weight": "distance"} if if_betweenness else {"distance": "distance"}
        partition_genes = {}
        t = tqdm(partition_genes_.items())
        for i, genes in t:
            if if_betweenness:
                kwargs["k"] = min(G.subgraph(genes).order(), k)
            t.set_description(
                f"Processing cluster {i}, size={G.subgraph(genes).order()}"
            )
            top_len = min(limit_anno_until, len(genes))
            top_gene_scores = dict(
                sorted(
                    compute_centrality(G.subgraph(genes), **kwargs).items(),
                    key=lambda x: x[1],
                    reverse=True,
                )[:top_len]
            )
            # Renormalizing centrality scores between 1 and 100, and rounding them to use later when
            # displaying wordclouds (higher score - higher "frequency" or word size)
            norm_top_gene_scores = dict(
                zip(
                    top_gene_scores.keys(),
                    scale_int(list(top_gene_scores.values()), 1, 100),
                )
            )
            partition_genes[i] = norm_top_gene_scores
        print("Filtered genes for generating the function word cloud..")

    elif isinstance(central_genes, dict):
        partition_genes = central_genes

    else:
        partition_genes = {
            i: {gene_: 1 for gene_ in gene_list}
            for i, gene_list in partition_genes_.items()
        }

    # If display gene function in the word clouds
    if display_func:
        # Computing functional annotation for each cluster as a concatenated list of annotations
        # Each annotation is weighted by its duplication gene_score times (e.g. a gene has score = 2 
        # -> the functional annotation is duplicated and have bigger font in WordCloud)
        partition_funcs = {
            i: " ".join(
                chain.from_iterable(
                    [
                        gene_func[gene_func.index == gene].to_list() * gene_score
                        for gene, gene_score in gene_score_list.items()
                    ]
                )
            )
            for i, gene_score_list in partition_genes.items()
        }

        # Generating word counts from aggregated gene annotation texts -> 
        # obtaining main (most frequent) function tokens
        word_counts = {
            i: WordCloud(
                max_words=30, min_font_size=15, stopwords=_STOPWORDS
            ).process_text(text)
            for i, text in partition_funcs.items()
        }
        word_counts = {
            i: (freqs if freqs else {"no found function": 1})
            for i, freqs in word_counts.items()
        }  # dealing with no word case
        wordclouds = {
            i: WordCloud(
                max_words=30,
                min_font_size=15,
                stopwords=_STOPWORDS,
                background_color="white",
                mask=get_elipsis_mask(),
            ).generate_from_frequencies(freqs)
            for i, freqs in word_counts.items()
        }

    # Display main genes in decreasing order of importance (top `top_len` genes)
    else:
        wordclouds = {
            i: WordCloud(
                max_words=30,
                min_font_size=15,
                background_color="white",
                mask=get_elipsis_mask(),
            ).generate_from_frequencies(gene_score_dict)
            for i, gene_score_dict in partition_genes.items()
        }

    # Plotting the word cloud
    partition_coords = {}
    for gene, coords in squeezed_pos.items():
        if partition[gene] not in partition_coords:
            partition_coords[partition[gene]] = [coords]
        else:
            partition_coords[partition[gene]] += [coords]
    for i, coords in partition_coords.items():
        x, y = zip(*coords)
        min_x, max_x = min(x), max(x)
        min_y, max_y = min(y), max(y)
        ax.imshow(
            wordclouds[i], interpolation="bilinear", extent=[min_x, max_x, min_y, max_y]
        )

    return ax
