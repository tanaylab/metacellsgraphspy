"""
Bar graphs of a metacells repository. See the Julia
`documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/bar_graphs.html>`__ for details.
"""

from typing import Optional

from dafpy import DafReader
from somegraphspy import SeriesBarsGraph

from .julia_import import _given
from .julia_import import jl

__all__ = [
    "degraded_genes_graph",
    "improved_genes_graph",
]


def improved_genes_graph(
    *,
    daf: DafReader,
    base_daf: DafReader,
    genes_count: Optional[int] = None,
    regulators_count: Optional[int] = None,
) -> SeriesBarsGraph:
    """
    The genes whose correlation with their cells the metacells improved in the most of the base neighborhoods. See the
    Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/bar_graphs.html#MetacellsGraphs.BarGraphs.improved_genes_graph>`__
    for details.
    """
    return SeriesBarsGraph.wrap_jl_object(
        jl.MetacellsGraphs.improved_genes_graph(
            daf=daf, base_daf=base_daf, **_given(genes_count=genes_count, regulators_count=regulators_count)
        )
    )


def degraded_genes_graph(
    *,
    daf: DafReader,
    base_daf: DafReader,
    genes_count: Optional[int] = None,
    regulators_count: Optional[int] = None,
) -> SeriesBarsGraph:
    """
    The genes whose correlation with their cells the metacells degraded in the most of the base neighborhoods. See the
    Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/bar_graphs.html#MetacellsGraphs.BarGraphs.degraded_genes_graph>`__
    for details.
    """
    return SeriesBarsGraph.wrap_jl_object(
        jl.MetacellsGraphs.degraded_genes_graph(
            daf=daf, base_daf=base_daf, **_given(genes_count=genes_count, regulators_count=regulators_count)
        )
    )
