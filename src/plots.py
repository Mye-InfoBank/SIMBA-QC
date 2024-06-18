from shiny import App, module, reactive, render, ui
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict
import numpy as np


@module.ui
def plots_ui():
    return ui.div(
        ui.card(ui.card_header("Number of cells"), ui.output_text("n_cells")),
        ui.card(ui.card_header("Scatter plot"), ui.output_plot("plot_scatter")),
        ui.card(
            ui.card_header("Histograms"),
            ui.output_ui("coloring_histograms"),
            ui.output_plot("plot_histograms"),
        ),
    )


@module.server
def plots_server(
    input,
    output,
    session,
    _adata: reactive.Value[ad.AnnData],
    _pretty_names: reactive.Value[Dict[str, str]],
    _distributions: reactive.Value[Dict[str, Dict[str, float]]],
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

        ax.scatter(
            adata.obs[x_col].astype(float),
            adata.obs[y_col].astype(float),
            c=adata.obs[color_col].astype(float),
            s=dot_size,
            cmap="viridis",
        )

        ax.set_xlabel(pretty_names[x_col])
        ax.set_ylabel(pretty_names[y_col])

        if adata.obs[color_col].nunique() > 1:
            # Add color bar
            cbar = plt.colorbar(ax.collections[0], ax=ax)
            cbar.set_label(pretty_names[color_col])

        return fig

    @output
    @render.ui
    def coloring_histograms():
        adata = _adata.get()

        if adata is None:
            return

        categorical_obs = [None] + [
            column
            for column in adata.obs.select_dtypes(
                include=["object", "category"]
            ).columns
        ]

        return ui.input_select("histo_coloring", "Coloring", categorical_obs)

    @output
    @render.plot
    def plot_histograms():
        adata = _adata.get()
        pretty_names = _pretty_names.get()
        distributions = _distributions.get()
        coloring = input["histo_coloring"].get()

        if adata is None:
            return

        n_plots = len(pretty_names)
        n_cols = min(2, n_plots)
        n_rows = int(np.ceil(n_plots / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols)

        for i, (col, pretty_name) in enumerate(pretty_names.items()):
            ax = axes[i // n_cols, i % n_cols] if n_rows > 1 else axes[i]

            kwargs = {"x": col, "ax": ax, "bins": 50}
            if coloring:
                kwargs["hue"] = coloring

            sns.histplot(adata.obs, **kwargs)
            ax.set_xlabel(pretty_name)

            current_distribution = distributions[col]

            # Add horizontal line for median
            ax.axvline(current_distribution["median"], color="r", linestyle="--")

            for mads in [1, 2, 3]:
                ax.axvline(
                    min(
                        current_distribution["median"]
                        + mads * current_distribution["std"],
                        current_distribution["max"],
                    ),
                    color="g",
                    linestyle=":",
                )
                ax.axvline(
                    max(
                        current_distribution["median"]
                        - mads * current_distribution["std"],
                        current_distribution["min"],
                    ),
                    color="g",
                    linestyle=":",
                )

        return fig

    @output
    @render.text
    def n_cells():
        adata = _adata.get()

        if adata is not None:
            n_cells = int(adata.n_obs)
            print(type(n_cells))
            return f"{n_cells} cells pass the current filtering thresholds"
