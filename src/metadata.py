from shiny import App, module, reactive, render, ui
import anndata as ad
import pandas as pd

mandatory_columns = ["batch", "cell_type", "condition", "sex", "patient", "tissue"]

@module.ui
def metadata_ui():
    return ui.div(
        ui.output_ui("column_cards")
    )

@module.server
def metadata_server(input, output, session,
                   _adata: reactive.Value[ad.AnnData],
                   _metadata: reactive.Value[pd.DataFrame]
                   ):
    _additional_columns = reactive.value([])
    _all_columns = reactive.value(mandatory_columns)

    @reactive.effect
    @reactive.event((input["add_column"]))
    def add_column():
        column = input["column_name"].get()
        if column and column not in _all_columns.get():
            _additional_columns.set(_additional_columns.get() + [column])
            ui.update_text("column_name", value="")

    @reactive.effect
    def update_columns():
        _all_columns.set(mandatory_columns + _additional_columns.get())

    @render.ui
    def column_cards():
        add_card = ui.card(
            ui.card_header("Add column"),
            ui.input_text("column_name", "", placeholder="Column name"),
            ui.input_action_button("add_column", "Add column")
        )

        return ui.layout_columns(
            *[ui.card(
                    ui.card_header(column),
                    ui.p("This is a column")
                )
            for column in _all_columns.get()],
            add_card
        )
