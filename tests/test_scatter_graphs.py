"""
Test the scatter graphs.

The Julia side tests what the graphs contain; what is worth testing here is that the wrapper reaches it and hands back
something the notebook can display - which is a ``somegraphspy`` graph rather than a raw Julia object.
"""

# The base repository is built the same way as in ``test_bar_graphs.py``, since both graphs are read against the same
# shape of base. The same disable is in that file, since the check needs both sides to report it.
# pylint: disable=duplicate-code

import dafpy as dp
import numpy as np
from somegraphspy import PointsGraph

import metacellsgraphspy as mg


def _test_daf() -> dp.DafWriter:
    """
    A repository with the least a gene-gene graph needs, and no types.
    """
    daf = dp.memory_daf(name="test!")
    daf.add_axis("gene", ["A", "B"])
    daf.add_axis("metacell", ["M1", "M2", "M3"])
    daf.add_axis("block", ["B1", "B2"])
    # Column-major, which is what ``Daf`` stores: the layout is not a detail it papers over.
    daf.set_matrix(
        "gene", "metacell", "linear_fraction", np.array([[0.1, 0.2, 0.3], [0.3, 0.2, 0.1]], dtype="float32", order="F")
    )
    daf.set_matrix(
        "gene", "block", "linear_fraction", np.array([[0.15, 0.35], [0.35, 0.15]], dtype="float32", order="F")
    )
    daf.set_vector("metacell", "umap_x", np.array([0.0, 1.0, 2.0], dtype="float32"))
    daf.set_vector("metacell", "umap_y", np.array([2.0, 1.0, 0.0], dtype="float32"))
    daf.set_vector("block", "umap_x", np.array([0.5, 1.5], dtype="float32"))
    daf.set_vector("block", "umap_y", np.array([1.5, 0.5], dtype="float32"))
    return daf


def test_metacells_gene_gene_graph() -> None:
    """
    The metacells graph is a points graph with a point per metacell.
    """
    graph = mg.metacells_gene_gene_graph(_test_daf(), x_gene="A", y_gene="B")
    assert isinstance(graph, PointsGraph)
    assert list(graph.data.points_xs) == [0.1, 0.2, 0.3]
    assert graph.figure is not None


def test_blocks_gene_gene_graph() -> None:
    """
    The blocks graph is a points graph with a point per block.
    """
    graph = mg.blocks_gene_gene_graph(_test_daf(), x_gene="A", y_gene="B")
    assert isinstance(graph, PointsGraph)
    assert list(graph.data.points_xs) == [0.15, 0.35]


def test_metacells_umap_graph() -> None:
    """
    The metacells UMAP graph is a points graph with a point per metacell.
    """
    graph = mg.metacells_umap_graph(_test_daf())
    assert isinstance(graph, PointsGraph)
    assert list(graph.data.points_xs) == [0.0, 1.0, 2.0]
    assert list(graph.data.points_ys) == [2.0, 1.0, 0.0]
    assert graph.figure is not None


def test_blocks_umap_graph() -> None:
    """
    The blocks UMAP graph is a points graph with a point per block.
    """
    graph = mg.blocks_umap_graph(_test_daf())
    assert isinstance(graph, PointsGraph)
    assert list(graph.data.points_xs) == [0.5, 1.5]


def test_gene_base_delta_correlations_graph() -> None:
    """
    The gene change graph is a points graph with a point per correlated base block.
    """
    correlation = "correlation_between_base_neighborhood_cells_and_punctuated_metacells"

    base_daf = dp.memory_daf(name="base!")
    base_daf.add_axis("gene", ["A", "B"])
    base_daf.add_axis("block", ["B1", "B2"])
    base_daf.add_axis("base_block", ["B1", "B2"])
    # Column-major, which is what ``Daf`` stores: the layout is not a detail it papers over.
    base_daf.set_matrix(
        "gene", "block", "linear_fraction", np.array([[0.1, 0.2], [0.3, 0.4]], dtype="float32", order="F")
    )
    base_daf.set_matrix(
        "gene", "base_block", correlation, np.array([[0.3, 0.7], [0.0, 0.0]], dtype="float32", order="F")
    )

    daf = dp.memory_daf(name="test!")
    daf.add_axis("gene", ["A", "B"])
    daf.add_axis("base_block", ["B1", "B2"])
    daf.set_matrix("gene", "base_block", correlation, np.array([[0.4, 0.6], [0.0, 0.0]], dtype="float32", order="F"))

    graph = mg.gene_base_delta_correlations_graph(daf=daf, base_daf=base_daf, gene="A")
    assert isinstance(graph, PointsGraph)
    assert [round(x, 3) for x in graph.data.points_xs] == [0.1, -0.1]
    assert list(graph.data.points_ys) == [0.1, 0.2]
    assert graph.figure is not None


def test_gene_fraction_regularization() -> None:
    """
    The regularization reaches the axes it is applied on.
    """
    graph = mg.metacells_gene_gene_graph(_test_daf(), x_gene="A", y_gene="B", gene_fraction_regularization=1e-3)
    assert graph.configuration.x_axis.log_regularization == 1e-3
    assert graph.configuration.y_axis.log_regularization == 1e-3
