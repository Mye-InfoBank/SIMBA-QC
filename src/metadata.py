from shiny import module, reactive, render, ui
import anndata as ad
import pandas as pd

mandatory_columns = ["batch", "cell_type", "condition", "sex", "patient", "tissue"]

@module.ui
def metadata_ui():
    return ui.navset_tab(
        ui.nav_panel("Columns",
                           ui.output_ui("column_cards")),
        ui.nav_panel("Preview", ui.output_data_frame("metadata_table"))
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
        type_string = input[select_accession].get() if select_accession in input else "Constant value"
        adata = _adata.get()

        if adata is None:
            print("No adata")
            return None

        if type_string == "Concat existing columns":
            select_concat_accession = f"select_concat_{column}"
            included_columns = list(input[select_concat_accession].get()) if select_concat_accession in input else []
            interface = ui.input_selectize(select_concat_accession, "Columns", _input_columns.get(), selected=included_columns, multiple=True)
        elif type_string == "Map existing":
            mapcol_accession = f"select_mapcol_{column}"
            available_columns = _input_columns.get()
            colselect = ui.input_select(mapcol_accession, "Column", available_columns)
            
            series = adata.obs[input[mapcol_accession].get() if mapcol_accession in input else available_columns[0]]
            unique_values = series.unique() if series is not None else []

            interface = ui.div(
                    colselect,
                    ui.accordion(
                        ui.accordion_panel("Mapping",
                            *[ui.input_text(f"mapping_{column}_{value}", 
                                            value, 
                                            input[f"mapping_{column}_{value}"].get() if f"mapping_{column}_{value}" in input else "", 
                                            placeholder="Unknown") for value in unique_values]
                    ))
                )
        elif type_string == "Constant value":
            value_accession = f"select_constant_{column}"
            existing_value = input[value_accession].get() if value_accession in input else ""
            interface = ui.input_text(value_accession, "Value", existing_value, placeholder="Unknown")
        else:
            raise ValueError(f"Unknown type {type_string}")


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
    
    @reactive.effect
    def update_metadata():
        adata = _adata.get()
        if adata is None:
            return
        metadata = pd.DataFrame(index=adata.obs.index)
        for column in _all_columns.get():
            select_accession = f"select_type_{column}"
            type_string = input[select_accession].get() if select_accession in input else None
            if type_string == "Concat existing columns":
                included_columns = list(input[f"select_concat_{column}"].get()) if f"select_concat_{column}" in input else []
                metadata[column] = adata.obs[included_columns].astype(str).apply(lambda x: "_".join(x), axis=1) if included_columns else "Unknown"
            elif type_string == "Map existing":
                mapcol_accession = f"select_mapcol_{column}"
                available_columns = _input_columns.get()
                series = adata.obs[input[mapcol_accession].get() if mapcol_accession in input else available_columns[0]]
                
                mapping = {value: input[f"mapping_{column}_{value}"].get() for value in series.unique()}
                metadata[column] = series.map(mapping)
            else:
                constant_accession = f"select_constant_{column}"
                metadata[column] = input[constant_accession].get() if constant_accession in input else "Unknown"
        _metadata.set(metadata)

    @render.data_frame
    def metadata_table():
        return _metadata.get()