from shiny import App, reactive, ui, render
import scanpy as sc
import os
import tempfile
import anndata as ad

from sliders import slider_ui, slider_server
from distributions import distributions_server
from plots import plots_server, plots_ui
from helpers import calculate_qc_metrics

_adata: reactive.Value[ad.AnnData] = reactive.value(None)
_adata_filtered: reactive.Value[ad.AnnData] = reactive.value(None)

_file_name = reactive.value(None)

_pretty_names = reactive.value({
    'total_counts': 'Total counts',
    'n_genes_by_counts': 'Number of genes with counts',
    'pct_counts_mt': 'Percentage of counts from mitochondrial genes',
})

_distributions = reactive.value({})

app_ui = ui.page_sidebar(
    ui.sidebar(ui.div(
        ui.input_file("file_input", label="Upload your file", accept=".h5ad"),
        slider_ui("sliders"),
        ui.download_button("download", "Download filtered data")
    )),
    plots_ui("plots"),
    title="Quality control for SIMBA🦁",
    window_title="SIMBA🦁: QC"
)


def server(input, output, session):
    distributions_server("distributions", _adata, _pretty_names, _distributions)
    slider_server("sliders", _adata, _adata_filtered, _pretty_names, _distributions)
    plots_server("plots", _adata_filtered, _pretty_names, _distributions)

    @reactive.effect
    def load_adata():
        file = input["file_input"].get()
        if file is None:
            return
        if len(file) > 1:
            print("Only one file at a time")
            return

        used_file = file[0]
        _file_name.set(used_file["name"])
        adata = sc.read(used_file["datapath"])
        calculate_qc_metrics(adata)
        _adata.set(adata)

    @render.download(
        filename = lambda: _file_name.get().replace(".h5ad", "_filtered.h5ad"),
    )
    def download():
        adata_filtered = _adata_filtered.get()
        if adata_filtered is None:
            return
        
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as temp:
            adata_filtered.write(temp.name)
            return temp.name

app = App(app_ui, server, static_assets=os.path.join(os.path.dirname(__file__), "media"))
