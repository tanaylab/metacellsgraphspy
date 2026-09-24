"""
Scatter graphs of a metacells repository. See the Julia
`documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/scatter_graphs.html>`__ for details.
"""

from typing import Optional

from dafpy import DafReader
from somegraphspy import PointsGraph

from .data_sources import Entries
from .julia_import import _given
from .julia_import import _to_julia_array
from .julia_import import jl

__all__ = [
    "gene_base_delta_correlations_graph",
    "gene_gene_graph",
    "umap_graph",
]


def gene_gene_graph(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    x_gene: str,
    y_gene: str,
    entries: Optional[Entries] = None,
    gene_fraction_regularization: Optional[float] = None,
) -> PointsGraph:
    """
    The expression of ``x_gene`` against ``y_gene``, a point per each of the ``entries`` of the ``axis`` (by default,
    the metacells), on log scale. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/scatter_graphs.html#MetacellsGraphs.ScatterGraphs.gene_gene_graph>`__
    for details.
    """
    return PointsGraph.wrap_jl_object(
        jl.MetacellsGraphs.gene_gene_graph(
            daf,
            x_gene=x_gene,
            y_gene=y_gene,
            **_given(
                axis=axis,
                entries=_to_julia_array(entries),
                gene_fraction_regularization=gene_fraction_regularization,
            ),
        )
    )


def umap_graph(daf: DafReader, *, axis: Optional[str] = None, entries: Optional[Entries] = None) -> PointsGraph:
    """
    The 2D UMAP embedding, a point per each of the ``entries`` of the ``axis`` (by default, the metacells). See the
    Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/scatter_graphs.html#MetacellsGraphs.ScatterGraphs.umap_graph>`__
    for details.
    """
    return PointsGraph.wrap_jl_object(
        jl.MetacellsGraphs.umap_graph(daf, **_given(axis=axis, entries=_to_julia_array(entries)))
    )


def gene_base_delta_correlations_graph(
    *,
    daf: DafReader,
    base_daf: DafReader,
    gene: str,
    axis: Optional[str] = None,
    base_axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    gene_fraction_regularization: Optional[float] = None,
) -> PointsGraph:
    """
    What the metacells did to one gene, a point per each of the ``entries`` of the ``axis`` (by default, the base
    blocks): how much its correlation with the cells changed, against how much of the gene there is in the entry. See
    the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/scatter_graphs.html#MetacellsGraphs.ScatterGraphs.gene_base_delta_correlations_graph>`__
    for details.
    """
    return PointsGraph.wrap_jl_object(
        jl.MetacellsGraphs.gene_base_delta_correlations_graph(
            daf=daf,
            base_daf=base_daf,
            gene=gene,
            **_given(
                axis=axis,
                base_axis=base_axis,
                entries=_to_julia_array(entries),
                gene_fraction_regularization=gene_fraction_regularization,
            ),
        )
    )
