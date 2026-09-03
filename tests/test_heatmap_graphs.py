"""
Test the heatmap graphs.

The Julia side tests what the graphs contain; what is worth testing here is that the wrapper reaches it and hands back
something the notebook can display - which is a ``somegraphspy`` graph rather than a raw Julia object.
"""

import dafpy as dp
import numpy as np
from somegraphspy import HeatmapGraph

import metacellsgraphspy as mg


def _test_daf() -> dp.DafWriter:
    """
    A repository with the least a heatmap graph needs, and no types.
    """
    daf = dp.memory_daf(name="test!")
    daf.add_axis("gene", ["A", "B", "C"])
    daf.add_axis("metacell", ["M1", "M2"])
    daf.add_axis("block", ["B1", "B2"])
    daf.set_vector("metacell", "block", np.array(["B1", "B2"], dtype=str))
    daf.set_vector("gene", "marker_rank", np.array([1, 2, np.iinfo(np.uint32).max], dtype="uint32"))
    daf.set_vector("gene", "is_skeleton", np.array([True, False, False]))
    # Column-major, which is what ``Daf`` stores: the layout is not a detail it papers over.
    daf.set_matrix(
        "gene",
        "metacell",
        "log_linear_fraction",
        np.array([[1.0, 3.0], [3.0, 1.0], [0.0, 0.0]], dtype="float32", order="F"),
    )
    daf.set_matrix(
        "gene",
        "block",
        "log_linear_fraction",
        np.array([[1.0, 3.0], [3.0, 1.0], [0.0, 0.0]], dtype="float32", order="F"),
    )
    return daf


def test_markers_metacells_heatmap_graph() -> None:
    """
    The metacells markers graph is a heatmap with a column per metacell, showing the fold from the median.
    """
    graph = mg.markers_metacells_heatmap_graph(_test_daf())
    assert isinstance(graph, HeatmapGraph)
    rows_names = graph.data.rows_names
    entries_values = graph.data.entries_values
    assert rows_names is not None and entries_values is not None
    assert list(rows_names) == ["A", "B"]
    assert list(entries_values[0, :]) == [-1.0, 1.0]
    assert graph.figure is not None


def test_markers_count() -> None:
    """
    The count reaches the choice of which markers are shown.
    """
    graph = mg.markers_metacells_heatmap_graph(_test_daf(), markers_count=1)
    rows_names = graph.data.rows_names
    assert rows_names is not None
    assert list(rows_names) == ["A"]


def test_skeletons_metacells_heatmap_graph() -> None:
    """
    The metacells skeletons graph shows the skeleton genes.
    """
    graph = mg.skeletons_metacells_heatmap_graph(_test_daf())
    assert isinstance(graph, HeatmapGraph)
    rows_names = graph.data.rows_names
    assert rows_names is not None
    assert list(rows_names) == ["A"]


def test_markers_blocks_heatmap_graph() -> None:
    """
    The blocks markers graph is a heatmap with a column per block.
    """
    graph = mg.markers_blocks_heatmap_graph(_test_daf())
    assert isinstance(graph, HeatmapGraph)
    columns_hovers = graph.data.columns_hovers
    assert columns_hovers is not None
    assert list(columns_hovers) == ["block: B1", "block: B2"]


def test_skeletons_blocks_heatmap_graph() -> None:
    """
    The blocks skeletons graph shows the skeleton genes.
    """
    graph = mg.skeletons_blocks_heatmap_graph(_test_daf())
    assert isinstance(graph, HeatmapGraph)
    rows_names = graph.data.rows_names
    assert rows_names is not None
    assert list(rows_names) == ["A"]


def test_group_by_block() -> None:
    """
    The grouping flags reach the groups of the columns.
    """
    graph = mg.markers_metacells_heatmap_graph(_test_daf(), group_by_block=True)
    columns_groups = graph.data.columns_groups
    assert columns_groups is not None
    assert list(columns_groups) == ["B1", "B2"]
