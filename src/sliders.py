from shiny import module, reactive, render, ui
import anndata as ad
from typing import Dict
import numpy as np
from helpers import calculate_qc_metrics


@module.ui
def slider_ui():
    return ui.div(
        ui.output_ui("button"),
        ui.output_ui("slider_sample"),
        ui.output_ui("slider_filters"),
        ui.output_ui("filter_button"),
    )


@module.server
def slider_server(
    input,
    output,
    session,
    _adata_meta: reactive.Value[ad.AnnData],
    _adata_qc: reactive.Value[ad.AnnData],
    _adata_filtered: reactive.Value[ad.AnnData],
    _pretty_names: reactive.Value[Dict[str, str]],
    _distributions: reactive.Value[Dict[str, Dict[str, float]]],
    _calculate_metrics_bool: reactive.value[bool],
):
    _adata_sample = reactive.value(None)
    _prev_mads = reactive.value({})
    _cur_filters = reactive.value({})
    _prev_filters = reactive.value({})

    @output
    @render.ui
    def slider_sample():
        adata = _adata_qc.get()
        if adata is None:
            return

        n_obs = adata.n_obs
        return ui.input_slider(
            "random_sample_size",
            "Random sample size",
            min(100, n_obs),
            n_obs,
            min(10000, n_obs),
            post=" cells",
        )

    @output
    @render.ui
    def button():
        calculate_metrics_bool = _calculate_metrics_bool.get()
        if calculate_metrics_bool is None:
            return

        if calculate_metrics_bool:
            return ui.input_task_button(
                "calculate_button",
                "Recalculate QC metrics",
                style="background-color: rgb(153, 0, 255); border-color: rgb(153, 0, 255);",
            )
        else:
            return None

    @ui.bind_task_button(button_id="calculate_button")
    @reactive.extended_task
    async def recalc_logic(adata_meta):
        adata_meta_copy = adata_meta.copy()
        calculate_qc_metrics(adata_meta_copy)
        return adata_meta_copy

    @reactive.effect
    @reactive.event(input.calculate_button, ignore_none=True)
    def handle_click():
        adata_meta = _adata_meta.get()
        if adata_meta is None:
            return
        recalc_logic(adata_meta)

    @reactive.effect
    def return_of_adata_meta():
        result = recalc_logic.result()
        _adata_qc.set(result)
        _calculate_metrics_bool.set(False)

    @reactive.effect
    def random_sample():
        adata = _adata_qc.get()
        if adata is None:
            return

        sample_size = input["random_sample_size"].get()

        adata_sample = adata[
            np.random.choice(
                adata.obs.index, min(sample_size, len(adata.obs)), replace=False
            )
        ]
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
            if distributions[col]["min"] == distributions[col]["max"]:
                continue
            else:
                mads = ui.input_slider(f"{col}_mads", "MADs", 0.25, 10, 2, step=0.25)
                absolute = ui.input_slider(
                    f"{col}_absolute",
                    "Absolute value",
                    distributions[col]["min"],
                    distributions[col]["max"],
                    [distributions[col]["min"], distributions[col]["max"]],
                )
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
            if distributions[col]["min"] == distributions[col]["max"]:
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

                min_val = (
                    distributions[col]["median"] - mads * distributions[col]["std"]
                )
                max_val = (
                    distributions[col]["median"] + mads * distributions[col]["std"]
                )

                ui.update_slider(absolute_name, value=[min_val, max_val])

                prev_mads[mads_name] = mads

    @reactive.effect
    def update_filters():
        filters = {}
        pretty_names = _pretty_names.get()

        for col in pretty_names.keys():
            input_name = f"{col}_absolute"
            if not input_name in input:
                continue
            min_val, max_val = input[input_name].get()
            filters[col] = (min_val, max_val)

        _cur_filters.set(filters)

    @render.ui
    def filter_button():
        filters = _cur_filters.get()
        prev_filters = _prev_filters.get()

        if filters == prev_filters:
            return None
        else:
            return ui.input_task_button(
                "filter",
                "Filter",
                style="background-color: rgb(153, 0, 255); border-color: rgb(153, 0, 255);",
            )

    @reactive.effect
    @reactive.event(input["filter"])
    def on_filter():
        adata = _adata_sample.get()
        filters = _cur_filters.get()

        if adata is None:
            return

        filter(adata, filters)

    @ui.bind_task_button(button_id="filter")
    @reactive.extended_task
    async def filter(adata: ad.AnnData, filters: dict):
        def filter_row(row):
            for col, thresholds in filters.items():
                min_val, max_val = thresholds
                if not (min_val <= row[col] <= max_val):
                    return False
            return True

        adata_filtered = adata[adata.obs.apply(filter_row, axis=1)]
        return adata_filtered, filters

    @reactive.effect
    def update_filtered():
        adata_filtered, filters = filter.result()
        if adata_filtered is not None:
            _adata_filtered.set(adata_filtered)
            _prev_filters.set(filters)
