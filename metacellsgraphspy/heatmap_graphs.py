"""
Heatmap graphs of a metacells repository. See the Julia
`documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/heatmap_graphs.html>`__ for details.
"""

from typing import Optional

from dafpy import DafReader
from somegraphspy import HeatmapGraph

from .data_sources import Entries
from .julia_import import _given
from .julia_import import _to_julia_array
from .julia_import import jl

__all__ = [
    "genes_heatmap_graph",
]


def genes_heatmap_graph(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    genes: Optional[Entries] = None,
    entries: Optional[Entries] = None,
    rows_axis_title: Optional[str] = None,
    columns_axis_title: Optional[str] = None,
    max_fold: Optional[float] = None,
) -> HeatmapGraph:
    """
    The expression of the ``genes`` (by default, all of them) in each of the ``entries`` of the ``axis`` (by default,
    the metacells), as the fold factor of the gene from its median across the entries. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/heatmap_graphs.html#MetacellsGraphs.HeatmapGraphs.genes_heatmap_graph>`__
    for details.
    """
    return HeatmapGraph.wrap_jl_object(
        jl.MetacellsGraphs.genes_heatmap_graph(
            daf,
            **_given(
                axis=axis,
                genes=_to_julia_array(genes),
                entries=_to_julia_array(entries),
                rows_axis_title=rows_axis_title,
                columns_axis_title=columns_axis_title,
                max_fold=max_fold,
            ),
        )
    )
