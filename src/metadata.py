from shiny import module, reactive, render, ui
import anndata as ad
import pandas as pd

mandatory_columns = ["batch", "cell_type", "condition", "sex", "patient", "tissue"]

@module.ui
def metadata_ui():
    return ui.accordion(
        ui.accordion_panel("Columns",
                           ui.output_ui("column_cards")),
        ui.accordion_panel("Preview")
    )

@module.server
def metadata_server(input, output, session,
                   _adata: reactive.Value[ad.AnnData],
                   _metadata: reactive.Value[pd.DataFrame]
                   ):
    _input_columns = reactive.value([])
    _additional_columns = reactive.value([])
    _all_columns = reactive.value(mandatory_columns)

    @reactive.effect
    def update_input_columns():
        adata = _adata.get()
        if adata is None:
            return
        categorical_columns = adata.obs.select_dtypes(include=["category"]).columns.tolist()
        _input_columns.set(categorical_columns)

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


    def column_card(column: str):
        select_accession = f"select_type_{column}"
        type_string = input[select_accession].get() if select_accession in input else None
        adata = _adata.get()

        if adata is None:
            print("No adata")
            return None

        if type_string == "Concat existing columns":
            interface = ui.input_selectize(f"select_{column}", "Columns", _input_columns.get(), multiple=True)
        elif type_string == "Map existing":
            mapcol_accession = f"select_mapcol_{column}"
            colselect = ui.input_select(mapcol_accession, "Column", _input_columns.get())
            
            series = adata.obs[input[mapcol_accession].get()] if mapcol_accession in input else None

            interface = ui.div(
                    colselect,
                    ui.accordion(
                        ui.accordion_panel("Mapping",
                            *([ui.input_text(f"mapping_{column}_{col}", col, placeholder="New value") for col in series.unique()] if series is not None else [])
                    ))
                )
        else:
            interface = ui.input_text(f"input_{column}", "Value", placeholder="Unknown")


        return ui.card(
                    ui.card_header(column),
                    ui.input_select(select_accession, "Type", ["Constant value", "Concat existing columns", "Map existing"], selected=type_string),
                    interface,
                    ui.card_footer(
                        ui.input_action_button("remove_column", "Remove")
                    )
                )

    @render.ui
    def column_cards():
        add_card = ui.card(
            ui.card_header("Add column"),
            ui.input_text("column_name", "", placeholder="Column name"),
            ui.card_footer(
                ui.input_action_button("add_column", "Add")
            )
        )

        return ui.layout_columns(
            *[column_card(column) for column in _all_columns.get()],
            add_card
        )
    

