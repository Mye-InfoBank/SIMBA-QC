from shiny import App, module, reactive, render, ui
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict
import numpy as np

@module.ui
def plots_ui():
    return ui.div(
        ui.card(
            ui.card_header("Scatter plot"),
            ui.output_plot("plot_scatter")
        ),
        ui.card(
            ui.card_header("Violin plots"),
            ui.output_plot("plot_violines")
        )
    )


@module.server
def plots_server(input, output, session,
                   _adata: reactive.Value[ad.AnnData],
                   _pretty_names: reactive.Value[Dict[str, str]],
                   _distributions: reactive.Value[Dict[str, Dict[str, float]]]
                   ):
    @output
    @render.plot
    def plot_scatter():
        adata = _adata.get()
        pretty_names = _pretty_names.get()

        if adata is None:
            return

        x_col = "n_genes_by_counts"
        y_col = "total_counts"
        color_col = "pct_counts_mt"
        dot_size = 1

        fig, ax = plt.subplots()

        ax.scatter(adata.obs[x_col],
                   adata.obs[y_col],
                   c=adata.obs[color_col],
                   s=dot_size,
                   cmap='viridis')
        
        ax.set_xlabel(pretty_names[x_col])
        ax.set_ylabel(pretty_names[y_col])

        if adata.obs[color_col].nunique() > 1:
            # Add color bar
            cbar = plt.colorbar(ax.collections[0], ax=ax)
            cbar.set_label(pretty_names[color_col])

        return fig

    @output
    @render.plot
    def plot_violines():
        adata = _adata.get()
        pretty_names = _pretty_names.get()
        distributions = _distributions.get()

        if adata is None:
            return

        n_plots = len(pretty_names)
        n_cols = min(2, n_plots)
        n_rows = int(np.ceil(n_plots / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols)

        for i, (col, pretty_name) in enumerate(pretty_names.items()):
            ax = axes[i // n_cols, i % n_cols] if n_rows > 1 else axes[i]
            sns.violinplot(adata.obs[col], ax=ax)
            ax.set_xlabel(pretty_name)
            ax.set_ylabel('Density')

            current_distribution = distributions[col]

            # Add horizontal line for median
            ax.axhline(current_distribution['median'], color='r', linestyle='--')

            for mads in [1, 2, 3]:
                ax.axhline(min(current_distribution['median'] + mads * current_distribution['std'], current_distribution['max']), color='g', linestyle=':')
                ax.axhline(max(current_distribution['median'] - mads * current_distribution['std'], current_distribution['min']), color='g', linestyle=':')

        return fig