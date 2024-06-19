from shiny import App, reactive, ui, render, Session
import scanpy as sc
import os
import tempfile
import anndata as ad
import json
import numpy as np

from sliders import slider_ui, slider_server
from distributions import distributions_server
from plots import plots_server, plots_ui
from metadata import metadata_server, metadata_ui
from helpers import calculate_qc_metrics


def convert_float32_to_float(data):
    if isinstance(data, dict):
        #print('dict')
        return {k: convert_float32_to_float(v) for k, v in data.items()}
    elif isinstance(data, list):
        #print('list')
        return [convert_float32_to_float(v) for v in data]
    elif isinstance(data, np.float32) or isinstance(data, np.float64):
        #print(f"Converting {data} to float")
        return float(data)
    elif isinstance(data, np.ndarray):
        #print('ndarray')
        return data.tolist()
    return data

# Monkey-patch json.dumps
original_json_dumps = json.dumps

def custom_json_dumps(data, *args, **kwargs):
    non_float_32_data = convert_float32_to_float(data)
    return original_json_dumps(non_float_32_data, *args, **kwargs)

json.dumps = custom_json_dumps

_pretty_names = reactive.value({
    'total_counts': 'Total counts',
    'n_genes_by_counts': 'Number of genes with counts',
    'pct_counts_mt': 'Percentage of counts from mitochondrial genes',
})

app_ui = ui.page_navbar(
    ui.nav_panel("1. Upload", ui.input_file("file_input", label="Upload your file", accept=".h5ad")),
    ui.nav_panel("2. Metadata", metadata_ui("metadata")),
    ui.nav_panel("3. Quality control",
        ui.layout_sidebar(
            ui.sidebar(
                slider_ui("sliders")
                ),
            plots_ui("plots")
        )
    ),
    ui.nav_panel("4. Download",
                 ui.download_button("download", "Download filtered data")
                 ),
    title="Preprocessing for SIMBA🦁",
    window_title="SIMBA🦁: Preprocessing"
)


def server(input, output, session: Session):
    _adata: reactive.Value[ad.AnnData] = reactive.value(None)
    _adata_meta: reactive.Value[ad.AnnData] = reactive.value(None)
    _adata_qc: reactive.Value[ad.AnnData] = reactive.value(None)
    _adata_filtered: reactive.Value[ad.AnnData] = reactive.value(None)
    _file_name = reactive.value(None)
    _distributions = reactive.value({})
    _metadata = reactive.value(None)
    _calculate_metrics_bool = reactive.value(False)

    distributions_server("distributions", _adata_qc, _pretty_names, _distributions)
    slider_server("sliders", _adata_meta, _adata_qc, _adata_filtered, _pretty_names, _distributions, _calculate_metrics_bool)
    plots_server("plots", _adata_filtered, _pretty_names, _distributions)
    metadata_server("metadata", _adata, _metadata)

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
        adata = sc.read_h5ad(used_file["datapath"])
        _adata.set(adata)
        _calculate_metrics_bool.set(True)

    @reactive.effect
    def update_adata_meta():
        adata = _adata.get()
        metadata = _metadata.get()
        if adata is None or metadata is None:
            return
        adata_meta = adata.copy()
        adata_meta.obs = metadata
        _adata_meta.set(adata_meta)
        _calculate_metrics_bool.set(True)

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
