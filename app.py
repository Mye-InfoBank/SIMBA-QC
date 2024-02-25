from shiny import App, reactive, ui
from modules import slider_ui, slider_server, distributions_server, plots_server, plots_ui
import scanpy as sc
from modules.helpers import calculate_qc_metrics
import os

_adata = reactive.value(None)
_adata_filtered = reactive.value(None)

_pretty_names = reactive.value({
    'total_counts': 'Total counts',
    'n_genes_by_counts': 'Number of genes with counts',
    'pct_counts_mt': 'Percentage of counts from mitochondrial genes',
})

_distributions = reactive.value({})

app_ui = ui.page_sidebar(
    ui.sidebar(ui.div(
        ui.input_file("file_input", label="Upload your file", accept=".h5ad"),
        slider_ui("sliders")
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
        print("Loading adata")
        file = input["file_input"].get()
        if file is None:
            return
        if len(file) > 1:
            print("Only one file at a time")
            return

        print(f"Reading file: {file}")
        adata = sc.read(file[0]["datapath"])
        calculate_qc_metrics(adata)
        _adata.set(adata)
        print("Loaded adata")

app = App(app_ui, server, static_assets=os.path.join(os.path.dirname(__file__), "media"))
