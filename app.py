from shiny import App, reactive, ui
from modules import slider_ui, slider_server, distributions_server, plots_server, plots_ui
import scanpy as sc
from modules.helpers import calculate_qc_metrics
import os

adata = sc.read_h5ad('minimal_example.h5ad')
calculate_qc_metrics(adata)

_adata = reactive.value(adata)
_adata_filtered = reactive.value(adata)

_pretty_names = reactive.value({
    'total_counts': 'Total counts',
    'n_genes_by_counts': 'Number of genes with counts',
    'pct_counts_mt': 'Percentage of counts from mitochondrial genes',
})

_distributions = reactive.value({})

app_ui = ui.page_sidebar(
    ui.sidebar(slider_ui("sliders")),
    plots_ui("plots"),
    title="Quality control for SIMBA🦁",
    window_title="SIMBA🦁: QC"
)


def server(input, output, session):
    distributions_server("distributions", _adata, _pretty_names, _distributions)
    slider_server("sliders", _adata, _adata_filtered, _pretty_names, _distributions)
    plots_server("plots", _adata_filtered, _pretty_names, _distributions)

app = App(app_ui, server, static_assets=os.path.join(os.path.dirname(__file__), "media"))
