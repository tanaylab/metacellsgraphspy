"""
Graphs of a metacells ``AnnData`` (e.g., read from an ``h5ad`` file).

These are the ``AnnData`` counterparts of some of the ``Daf`` graphs. They are implemented in Python only, using the
:py:mod:`~metacellsgraphspy.anndata_data_sources`. Each ``obs`` entry is a metacell, each ``var`` entry is a gene, and
``X`` holds the expression level of each gene in each metacell (typically a geometric mean, used as it is stored).
"""

from typing import Optional

from anndata import AnnData  # type: ignore
from somegraphspy import HeatmapGraph
from somegraphspy import PointsGraph
from somegraphspy import heatmap_graph
from somegraphspy import points_graph

from .anndata_data_sources import ad_fill_gene_expression
from .anndata_data_sources import ad_fill_genes_fold_matrix
from .data_sources import Entries
from .heatmap_graphs import put_genes_heatmap_configuration

__all__ = [
    "ad_gene_gene_graph",
    "ad_genes_heatmap_graph",
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


def ad_genes_heatmap_graph(
    adata: AnnData,
    *,
    genes: Optional[Entries] = None,
    entries: Optional[Entries] = None,
    rows_axis_title: Optional[str] = None,
    columns_axis_title: Optional[str] = None,
    max_fold: Optional[float] = None,
) -> HeatmapGraph:
    """
    The expression of the ``genes`` (by default, all of them) in each of the ``entries`` of the ``adata`` (by default,
    all the metacells), as the fold factor of the gene from its median across the entries. For a markers heatmap, pick
    the genes with :py:func:`~metacellsgraphspy.anndata_data_sources.ad_get_top_marker_gene_indices`. This is the
    ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.heatmap_graphs.genes_heatmap_graph`.
    """
    graph = heatmap_graph()
    ad_fill_genes_fold_matrix(graph.entries_matrix_fields(), adata, genes=genes, entries=entries, max_fold=max_fold)
    matrix = graph.data.entries.matrix
    assert matrix is not None
    n_genes, n_entries = matrix.shape
    put_genes_heatmap_configuration(
        graph.configuration,
        n_genes=n_genes,
        n_entries=n_entries,
        rows_axis_title=rows_axis_title,
        columns_axis_title=columns_axis_title,
    )
    return graph
