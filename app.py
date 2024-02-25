from shiny import reactive
from shiny.express import input, render, ui
import scanpy as sc
from helpers import calculate_qc_metrics, clean
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

adata = sc.read_h5ad('minimal_example.h5ad')
calculate_qc_metrics(adata)

pretty_names = {
    'total_counts': 'Total counts',
    'n_genes_by_counts': 'Number of genes with counts',
    'pct_counts_mt': 'Max mito percentage',
}

def get_distribution(column: pd.Series):
    return {
        "min": column.min(),
        "max": column.max(),
        "median": column.median(),
        "std": column.std(),
    }

distributions = {
    col: get_distribution(adata.obs[col]) for col in pretty_names.keys()
}

ui.h3('Quality Control for SIMBA🦁')

with ui.sidebar():
    ui.input_slider('n_cells', 'Number of cells', min(100, adata.n_obs), adata.n_obs, min(10000, adata.n_obs))

    for col, pretty_name in pretty_names.items():
        if distributions[col]['min'] == distributions[col]['max']:
            continue
        else:
            ui.input_slider(col,
                            pretty_name,
                            distributions[col]['min'],
                            distributions[col]['max'],
                            [distributions[col]['min'],
                            distributions[col]['max']])
    
    ui.input_action_button('plot', 'Update plots')

_subset = reactive.value(adata)
_filtered = reactive.value(adata)

@reactive.effect
def subset():
    n_cells = input.n_cells.get()
    # Randomly sample n_cells cells
    adata_subset = adata[np.random.choice(adata.obs.index, n_cells), :]
    _subset.set(adata_subset)

@reactive.effect
def filter():
    def filter_row(row):
        for col in pretty_names.keys():
            if not col in input:
                continue
            min_val, max_val = input[col].get()
            if not (min_val <= row[col] <= max_val):
                return False
        return True
    adata_subset = _subset.get()

    adata_filtered = adata_subset[adata_subset.obs.apply(filter_row, axis=1), :]

    _filtered.set(adata_filtered)

@render.plot
@reactive.event(input.plot)
def plot_scatter():
    adata_filtered = _filtered.get()

    fig, ax = plt.subplots()
    sc.pl.scatter(adata_filtered,
                  x='n_genes_by_counts',
                  y='total_counts',
                  color='pct_counts_mt',
                  ax=ax,
                  show=False)

    return clean(fig)

@render.plot
@reactive.event(input.plot)
def plot_violines():
    adata_filtered = _filtered.get()

    n_plots = len(pretty_names)
    n_cols = 2
    n_rows = n_plots // n_cols + n_plots % n_cols
    fig, axes = plt.subplots(n_rows, n_cols)

    for i, (col, pretty_name) in enumerate(pretty_names.items()):
        ax = axes[i // n_cols, i % n_cols] if n_rows > 1 else axes[i % n_cols]
        sns.violinplot(y=col, data=adata_filtered.obs, ax=ax)
        ax.set_title(pretty_name)

        # Add horizontal line for the median
        ax.axhline(distributions[col]['median'], color='r', linestyle='--')

        for mads in [1, 2, 3]:
            ax.axhline(min(
                distributions[col]['median'] + mads * distributions[col]['std'],
                distributions[col]['max']),
                       color='r', linestyle=':', alpha=0.5)
            ax.axhline(max(
                distributions[col]['median'] - mads * distributions[col]['std'],
                distributions[col]['min']),
                       color='r', linestyle=':', alpha=0.5)
    
    return fig