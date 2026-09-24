"""
Test the heatmap graphs.

The Julia side tests what the graphs contain; what is worth testing here is that the wrapper reaches it and hands back
something the notebook can display - which is a ``somegraphspy`` graph rather than a raw Julia object.
"""

from typing import Any
from typing import List

import dafpy as dp
import numpy as np
from somegraphspy import HeatmapGraph

import metacellsgraphspy as mg


def _list(values: Any) -> List[Any]:
    """
    A graph field as a list, once it is known to be filled.
    """
    assert values is not None
    return list(values)


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


def test_genes_heatmap_graph() -> None:
    """
    The genes heatmap is a heatmap with a row per gene and a column per metacell, showing the fold from the median.
    """
    graph = mg.genes_heatmap_graph(_test_daf())
    assert isinstance(graph, HeatmapGraph)
    assert _list(graph.data.rows.entities.names) == ["A", "B", "C"]
    assert _list(graph.data.columns.entities.names) == ["M1", "M2"]
    matrix = graph.data.entries.matrix
    assert matrix is not None
    assert list(matrix[0, :]) == [-1.0, 1.0]
    assert graph.figure is not None


def test_genes_selection() -> None:
    """
    The gene selection sources pick which genes are shown.
    """
    daf = _test_daf()
    graph = mg.genes_heatmap_graph(daf, genes=mg.get_top_marker_gene_indices(daf, markers_count=1))
    assert _list(graph.data.rows.entities.names) == ["A"]

    graph = mg.genes_heatmap_graph(daf, genes=mg.get_skeleton_gene_indices(daf))
    assert _list(graph.data.rows.entities.names) == ["A"]

    graph = mg.genes_heatmap_graph(daf, genes=["B", "C"])
    assert _list(graph.data.rows.entities.names) == ["B", "C"]


def test_genes_heatmap_graph_of_blocks() -> None:
    """
    The axis picks which entries the columns are.
    """
    graph = mg.genes_heatmap_graph(_test_daf(), axis="block", columns_axis_title="Blocks!")
    assert isinstance(graph, HeatmapGraph)
    assert _list(graph.data.columns.entities.names) == ["B1", "B2"]
    assert graph.configuration.columns.title == "Blocks!"
