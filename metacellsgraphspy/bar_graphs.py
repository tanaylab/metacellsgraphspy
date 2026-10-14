"""
Bar graphs of a metacells repository. See the Julia
`documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/bar_graphs.html>`__ for details.
"""

from typing import Optional

from dafpy import DafReader
from somegraphspy import SeriesBarsGraph
from somegraphspy import SeriesBarsGraphConfiguration
from somegraphspy.julia_import import _to_julia

from .julia_import import _given
from .julia_import import jl

__all__ = [
    "add_changed_genes_default_annotations",
    "degraded_genes_graph",
    "improved_genes_graph",
    "put_changed_genes_configuration",
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


def put_changed_genes_configuration(configuration: SeriesBarsGraphConfiguration) -> None:
    """
    Lay out a butterfly of genes, as :py:func:`improved_genes_graph` and :py:func:`degraded_genes_graph` show it. See
    the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/bar_graphs.html#MetacellsGraphs.BarGraphs.put_changed_genes_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_changed_genes_configuration_b(_to_julia(configuration))


def add_changed_genes_default_annotations(graph: SeriesBarsGraph, *, daf: DafReader, base_daf: DafReader) -> None:
    """
    Add the default annotations of the bars (the genes) of an :py:func:`improved_genes_graph` or a
    :py:func:`degraded_genes_graph`, read from the ``base_daf``, leaving out those whose data it doesn't have. See the
    Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/bar_graphs.html#MetacellsGraphs.BarGraphs.add_changed_genes_default_annotations!>`__
    for details.
    """
    jl.MetacellsGraphs.add_changed_genes_default_annotations_b(_to_julia(graph), daf=daf, base_daf=base_daf)
