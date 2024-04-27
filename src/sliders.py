from shiny import App, module, reactive, render, ui
import anndata as ad
from typing import Dict
import numpy as np
import pandas as pd
from helpers import calculate_qc_metrics

@module.ui
def slider_ui():
    return ui.div(
        ui.output_ui("slider_sample"),
        ui.output_ui("slider_filters"),
        ui.input_task_button("calculate_button", "Calculate Metadata new")
    )


@module.server
def slider_server(input, output, session,
                   _adata: reactive.Value[ad.AnnData],
                   _metadata: reactive.value[pd.DataFrame],
                   _adata_filtered: reactive.Value[ad.AnnData],
                   _pretty_names: reactive.Value[Dict[str, str]],
                   _distributions: reactive.Value[Dict[str, Dict[str, float]]],
                   _calculate_metrics_bool: reactive.value[bool],
                   ):
    _adata_sample = reactive.value(None)
    _prev_mads = reactive.value({})
    
    @output
    @render.ui
    def slider_sample():
        adata = _adata.get()
        if adata is None:
            return

        n_obs = adata.n_obs
        return ui.input_slider('random_sample_size', 'Random sample size', min(100, n_obs), n_obs, min(10000, n_obs), post=" cells")
    
    @reactive.effect
    def calculate_qc_metrics_and_update():
        adata = _adata.get()
        adata_meta = adata.copy()
        metadata = _metadata.get()
        adata_meta.obs = metadata.copy()
        calculate_qc_metrics(adata_meta)
        _calculate_metrics_bool.set(False)
        #_adata.set(adata)
        #_adata_meta.set(adata_meta)
        
    @reactive.effect
    @reactive.event(input.calculate_button, ignore_none=False)
    def handle_click():
        calculate_qc_metrics_and_update()

    @reactive.effect
    def showing_button():
        if _calculate_metrics_bool.get():
            ui.show("calculate_button")
        else:
            ui.hide("calculate_button")

    @reactive.effect
    def random_sample():
        adata = _adata.get()
        sample_size = input['random_sample_size'].get()
        
        if adata is None:
            return

        adata_sample = adata[np.random.choice(adata.obs.index, sample_size, replace=False)]
        _adata_sample.set(adata_sample)

    @output
    @render.ui
    def slider_filters():
        pretty_names = _pretty_names.get()
        distributions = _distributions.get()

        if distributions == {}:
            return ui.div(ui.p("No data available"))

        panels = []

        for col, pretty_name in pretty_names.items():
            if distributions[col]['min'] == distributions[col]['max']:
                continue
            else:
                mads = ui.input_slider(f"{col}_mads", "MADs", 0.25, 10, 2, step=0.25)
                absolute = ui.input_slider(f"{col}_absolute",
                                "Absolute value",
                                distributions[col]['min'],
                                distributions[col]['max'],
                                [distributions[col]['min'],
                                distributions[col]['max']])
                panel = ui.accordion_panel(pretty_name, mads, absolute)
                panels.append(panel)        
        return ui.accordion(*panels)
    
    @reactive.effect
    def update_absolute_by_mads():
        pretty_names = _pretty_names.get()
        distributions = _distributions.get()
        prev_mads = _prev_mads.get()

        if distributions == {}:
            return

        for col, _ in pretty_names.items():
            if distributions[col]['min'] == distributions[col]['max']:
                continue
            else:
                mads_name = f"{col}_mads"
                absolute_name = f"{col}_absolute"

                if not mads_name in input or not absolute_name in input:
                    continue

                mads = input[mads_name].get()
                if mads_name in prev_mads and prev_mads[mads_name] == mads:
                    continue

                absolute = input[absolute_name].get()
                if mads is None or absolute is None:
                    continue

                min_val = distributions[col]['median'] - mads * distributions[col]['std']
                max_val = distributions[col]['median'] + mads * distributions[col]['std']

                ui.update_slider(absolute_name, value=[min_val, max_val])

                prev_mads[mads_name] = mads
    
    @reactive.effect
    def filter():
        def filter_row(row):
            for col in pretty_names.keys():
                input_name = f"{col}_absolute"
                if not input_name in input:
                    continue
                min_val, max_val = input[input_name].get()
                if not (min_val <= row[col] <= max_val):
                    return False
            return True

        adata = _adata_sample.get()
        pretty_names = _pretty_names.get()

        if adata is None:
            return

        adata_filtered = adata[adata.obs.apply(filter_row, axis=1)]
        _adata_filtered.set(adata_filtered)
