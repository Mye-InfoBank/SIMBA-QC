import anndata as ad
import scanpy as sc
from typing import List

def __annotate_genes__(adata: ad.AnnData, use_col: str | None = None):
    col = adata.var_names if use_col is None else adata.var[use_col]

    adata.var["mt"] = col.str.startswith("MT-")
    adata.var["ribo"] = col.str.startswith(("RPL", "RPS"))
    adata.var["hb"] = col.str.contains(("^HB[^(P)]"))


def calculate_qc_metrics(adata: ad.AnnData, percent_top: List[int] = [20], use_col: str | None = None):
    __annotate_genes__(adata, use_col=use_col)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        log1p=True,
        percent_top=percent_top,
    )

def clean(fig):
    while len(fig.get_axes()) > 1:
        fig.delaxes(fig.get_axes()[1])
    return fig