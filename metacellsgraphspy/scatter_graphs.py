"""
Scatter graphs of a metacells repository. See the Julia
`documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/scatter_graphs.html>`__ for details.
"""

from typing import Optional

from dafpy import DafReader
from somegraphspy import PointsGraph

from .julia_import import _given
from .julia_import import jl

__all__ = [
    "blocks_gene_gene_graph",
    "metacells_gene_gene_graph",
]


def metacells_gene_gene_graph(
    daf: DafReader,
    *,
    x_gene: str,
    y_gene: str,
    gene_fraction_regularization: Optional[float] = None,
) -> PointsGraph:
    """
    The expression of one gene against another, a point per metacell, on log scale. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/scatter_graphs.html#MetacellsGraphs.ScatterGraphs.metacells_gene_gene_graph>`__
    for details.
    """
    return PointsGraph.wrap_jl_object(
        jl.MetacellsGraphs.metacells_gene_gene_graph(
            daf,
            x_gene=x_gene,
            y_gene=y_gene,
            **_given(gene_fraction_regularization=gene_fraction_regularization),
        )
    )


def blocks_gene_gene_graph(
    daf: DafReader,
    *,
    x_gene: str,
    y_gene: str,
    gene_fraction_regularization: Optional[float] = None,
) -> PointsGraph:
    """
    The expression of one gene against another, a point per block, on log scale. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/scatter_graphs.html#MetacellsGraphs.ScatterGraphs.blocks_gene_gene_graph>`__
    for details.
    """
    return PointsGraph.wrap_jl_object(
        jl.MetacellsGraphs.blocks_gene_gene_graph(
            daf,
            x_gene=x_gene,
            y_gene=y_gene,
            **_given(gene_fraction_regularization=gene_fraction_regularization),
        )
    )
