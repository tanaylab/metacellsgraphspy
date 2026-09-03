"""
Heatmap graphs of a metacells repository. See the Julia
`documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/heatmap_graphs.html>`__ for details.
"""

from typing import Optional

from dafpy import DafReader
from somegraphspy import HeatmapGraph

from .julia_import import _given
from .julia_import import jl

__all__ = [
    "markers_blocks_heatmap_graph",
    "markers_metacells_heatmap_graph",
    "skeletons_blocks_heatmap_graph",
    "skeletons_metacells_heatmap_graph",
]


def markers_metacells_heatmap_graph(
    daf: DafReader,
    *,
    markers_count: Optional[int] = None,
    group_by_type: Optional[bool] = None,
    group_by_block: Optional[bool] = None,
    use_global_flow_order: Optional[bool] = None,
) -> HeatmapGraph:
    """
    The expression of the best marker genes in each metacell, as the fold factor of the gene from its median. See the
    Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/heatmap_graphs.html#MetacellsGraphs.HeatmapGraphs.markers_metacells_heatmap_graph>`__
    for details.
    """
    return HeatmapGraph.wrap_jl_object(
        jl.MetacellsGraphs.markers_metacells_heatmap_graph(
            daf,
            **_given(
                markers_count=markers_count,
                group_by_type=group_by_type,
                group_by_block=group_by_block,
                use_global_flow_order=use_global_flow_order,
            ),
        )
    )


def skeletons_metacells_heatmap_graph(
    daf: DafReader,
    *,
    group_by_type: Optional[bool] = None,
    group_by_block: Optional[bool] = None,
    use_global_flow_order: Optional[bool] = None,
) -> HeatmapGraph:
    """
    The expression of the skeleton genes in each metacell, as the fold factor of the gene from its median. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/heatmap_graphs.html#MetacellsGraphs.HeatmapGraphs.skeletons_metacells_heatmap_graph>`__
    for details.
    """
    return HeatmapGraph.wrap_jl_object(
        jl.MetacellsGraphs.skeletons_metacells_heatmap_graph(
            daf,
            **_given(
                group_by_type=group_by_type,
                group_by_block=group_by_block,
                use_global_flow_order=use_global_flow_order,
            ),
        )
    )


def markers_blocks_heatmap_graph(
    daf: DafReader,
    *,
    markers_count: Optional[int] = None,
    group_by_type: Optional[bool] = None,
    use_global_flow_order: Optional[bool] = None,
) -> HeatmapGraph:
    """
    The expression of the best marker genes in each block, as the fold factor of the gene from its median. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/heatmap_graphs.html#MetacellsGraphs.HeatmapGraphs.markers_blocks_heatmap_graph>`__
    for details.
    """
    return HeatmapGraph.wrap_jl_object(
        jl.MetacellsGraphs.markers_blocks_heatmap_graph(
            daf,
            **_given(
                markers_count=markers_count,
                group_by_type=group_by_type,
                use_global_flow_order=use_global_flow_order,
            ),
        )
    )


def skeletons_blocks_heatmap_graph(
    daf: DafReader,
    *,
    group_by_type: Optional[bool] = None,
    use_global_flow_order: Optional[bool] = None,
) -> HeatmapGraph:
    """
    The expression of the skeleton genes in each block, as the fold factor of the gene from its median. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/heatmap_graphs.html#MetacellsGraphs.HeatmapGraphs.skeletons_blocks_heatmap_graph>`__
    for details.
    """
    return HeatmapGraph.wrap_jl_object(
        jl.MetacellsGraphs.skeletons_blocks_heatmap_graph(
            daf,
            **_given(
                group_by_type=group_by_type,
                use_global_flow_order=use_global_flow_order,
            ),
        )
    )
