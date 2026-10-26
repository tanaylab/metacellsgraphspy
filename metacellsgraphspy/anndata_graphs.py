"""
Graphs of a metacells ``AnnData`` (e.g., read from an ``h5ad`` file).

These are the ``AnnData`` counterparts of some of the ``Daf`` graphs. They are implemented in Python only, using the
:py:mod:`~metacellsgraphspy.anndata_data_sources`. Each ``obs`` entry is a metacell, each ``var`` entry is a gene, and
``X`` holds the expression level of each gene in each metacell (typically a geometric mean, used as it is stored).
"""

from anndata import AnnData  # type: ignore
from somegraphspy import HeatmapGraph
from somegraphspy import PointsGraph
from somegraphspy import heatmap_graph
from somegraphspy import points_graph

from .anndata_data_sources import ad_fill_gene_expression
from .anndata_data_sources import ad_fill_genes_fold_matrix
from .anndata_data_sources import ad_fill_is_lateral
from .anndata_data_sources import ad_fill_type
from .anndata_data_sources import ad_get_top_marker_gene_indices
from .data_sources import TOP_MARKERS_COUNT_FOR_GRAPHS
from .data_sources import Entries
from .heatmap_graphs import put_genes_heatmap_configuration
from .julia_import import DEFAULT
from .julia_import import DefaultValue

__all__ = [
    "ad_add_genes_heatmap_default_columns_annotations",
    "ad_add_genes_heatmap_default_rows_annotations",
    "ad_gene_gene_graph",
    "ad_genes_heatmap_graph",
]


def ad_gene_gene_graph(
    adata: AnnData,
    *,
    x_gene: str,
    y_gene: str,
    entries: Entries | None | DefaultValue = DEFAULT,
    gene_fraction_regularization: float | DefaultValue = DEFAULT,
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
    genes: Entries | None | DefaultValue = DEFAULT,
    entries: Entries | None | DefaultValue = DEFAULT,
    rows_axis_title: str | None | DefaultValue = DEFAULT,
    columns_axis_title: str | None | DefaultValue = DEFAULT,
    max_fold: float | DefaultValue = DEFAULT,
) -> HeatmapGraph:
    """
    The expression of the ``genes`` (by default, the
    :py:obj:`~metacellsgraphspy.data_sources.TOP_MARKERS_COUNT_FOR_GRAPHS` strongest markers, see
    :py:func:`~metacellsgraphspy.anndata_data_sources.ad_get_top_marker_gene_indices`) in each of the ``entries`` of the
    ``adata`` (by default, all the metacells), as the fold factor of the gene from its median across the entries. This
    is the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.heatmap_graphs.genes_heatmap_graph`.
    """
    if genes is None or isinstance(genes, DefaultValue):
        genes = ad_get_top_marker_gene_indices(adata, markers_count=TOP_MARKERS_COUNT_FOR_GRAPHS)
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


def ad_add_genes_heatmap_default_rows_annotations(graph: HeatmapGraph, adata: AnnData) -> None:
    """
    Add the default annotations of the rows (the genes) of an :py:func:`ad_genes_heatmap_graph` of the ``adata``:
    whether each gene is lateral, if the ``var`` has a ``lateral_gene`` column. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.heatmap_graphs.add_genes_heatmap_default_rows_annotations`.
    """
    if "lateral_gene" in adata.var:
        ad_fill_is_lateral(graph.rows_annotations_colors_vector_fields(graph.add_rows_annotation()), adata)


def ad_add_genes_heatmap_default_columns_annotations(
    graph: HeatmapGraph, adata: AnnData, *, type_property: str, type_colors_csv: str
) -> None:
    """
    Add the default annotations of the columns (the metacells) of an :py:func:`ad_genes_heatmap_graph` of the
    ``adata``: the type of each metacell, in the ``type_property`` of the ``obs`` (if it has one), colored by the
    ``type_colors_csv`` file. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.heatmap_graphs.add_genes_heatmap_default_columns_annotations`.
    """
    if type_property in adata.obs:
        ad_fill_type(
            graph.columns_annotations_colors_vector_fields(graph.add_columns_annotation()),
            adata,
            type_property=type_property,
            type_colors_csv=type_colors_csv,
        )
