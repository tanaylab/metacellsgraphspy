"""
Graphs of a metacells ``AnnData`` (e.g., read from an ``h5ad`` file).

These are the ``AnnData`` counterparts of some of the ``Daf`` graphs. They are implemented in Python only, using the
:py:mod:`~metacellsgraphspy.anndata_data_sources`. Each ``obs`` entry is a metacell, each ``var`` entry is a gene, and
``X`` holds the fraction of the UMIs of each metacell in each gene.
"""

from typing import Optional

from anndata import AnnData  # type: ignore
from somegraphspy import PointsGraph
from somegraphspy import points_graph

from .anndata_data_sources import ad_fill_gene_expression
from .data_sources import Entries

__all__ = [
    "ad_gene_gene_graph",
]


def ad_gene_gene_graph(
    adata: AnnData,
    *,
    x_gene: str,
    y_gene: str,
    entries: Optional[Entries] = None,
    gene_fraction_regularization: Optional[float] = None,
) -> PointsGraph:
    """
    The expression of ``x_gene`` against ``y_gene``, a point per each of the ``entries`` of the ``adata``, on log scale.
    This is the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.scatter_graphs.gene_gene_graph`.
    """
    graph = points_graph()
    ad_fill_gene_expression(
        graph.x_axis_vector_fields(),
        adata,
        gene=x_gene,
        entries=entries,
        gene_fraction_regularization=gene_fraction_regularization,
    )
    ad_fill_gene_expression(
        graph.y_axis_vector_fields(),
        adata,
        gene=y_gene,
        entries=entries,
        gene_fraction_regularization=gene_fraction_regularization,
    )
    return graph
