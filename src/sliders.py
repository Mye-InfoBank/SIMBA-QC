from shiny import App, module, reactive, render, ui
import anndata as ad
from typing import Dict
import numpy as np

@module.ui
def slider_ui():
    return ui.div(
        ui.output_ui("slider_sample"),
        ui.output_ui("slider_filters")
    )


@module.server
def slider_server(input, output, session,
                   _adata: reactive.Value[ad.AnnData],
                   _adata_filtered: reactive.Value[ad.AnnData],
                   _pretty_names: reactive.Value[Dict[str, str]],
                   _distributions: reactive.Value[Dict[str, Dict[str, float]]]
                   ):
    _adata_sample = reactive.value(None)

    @output
    @render.text
    def out():
        return f"Click count is {_adata.get()}\n\nPretty names are {_pretty_names.get()}\n\nDistributions are {_distributions.get()}"
    
    @output
    @render.ui
    def slider_sample():
        adata = _adata.get()
        if adata is None:
            return

        n_obs = adata.n_obs
        return ui.input_slider('random_sample_size', 'Random sample size', min(100, n_obs), n_obs, min(10000, n_obs))
    
    @reactive.effect
    def random_sample():
        adata = _adata.get()
        if adata is None:
            return

        sample_size = input['random_sample_size'].get() if 'random_sample_size' in input else adata.n_obs

        adata_sample = adata[np.random.choice(adata.obs.index, sample_size, replace=False)]
        _adata_sample.set(adata_sample)

    @output
    @render.ui
    def slider_filters():
        pretty_names = _pretty_names.get()
        distributions = _distributions.get()

        if distributions == {}:
            return ui.div(ui.p("No data available"))

        sliders = []

        for col, pretty_name in pretty_names.items():
            if distributions[col]['min'] == distributions[col]['max']:
                continue
            else:
                slider = ui.input_slider(col,
                                pretty_name,
                                distributions[col]['min'],
                                distributions[col]['max'],
                                [distributions[col]['min'],
                                distributions[col]['max']])
                sliders.append(slider)
        
        return ui.div(*sliders)
    
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

        adata = _adata_sample.get()
        pretty_names = _pretty_names.get()

        if adata is None:
            return

        adata_filtered = adata[adata.obs.apply(filter_row, axis=1)]
        _adata_filtered.set(adata_filtered)