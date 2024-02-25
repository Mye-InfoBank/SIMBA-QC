from shiny import module, reactive
import anndata as ad
from typing import Dict
import pandas as pd

def get_distribution(column: pd.Series):
    return {
        "min": column.min(),
        "max": column.max(),
        "median": column.median(),
        "std": column.std(),
    }

@module.server
def distributions_server(input, output, session,
                   _adata: reactive.Value[ad.AnnData],
                   _pretty_names: reactive.Value[Dict[str, str]],
                   _distributions: reactive.Value[Dict[str, Dict[str, float]]]
                   ):

    @reactive.effect
    def update_distributions():
        adata = _adata.get()
        distributions = {
            col: get_distribution(adata.obs[col]) for col in _pretty_names.get().keys()
        }
        _distributions.set(distributions)