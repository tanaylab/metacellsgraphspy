"""
Test the scatter graphs.

The Julia side tests what the graphs contain; what is worth testing here is that the wrapper reaches it and hands back
something the notebook can display - which is a ``somegraphspy`` graph rather than a raw Julia object.
"""

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
    # Column-major, which is what `Daf` stores: the layout is not a detail it papers over.
    daf.set_matrix(
        "gene", "metacell", "linear_fraction", np.array([[0.1, 0.2, 0.3], [0.3, 0.2, 0.1]], dtype="float32", order="F")
    )
    daf.set_matrix(
        "gene", "block", "linear_fraction", np.array([[0.15, 0.35], [0.35, 0.15]], dtype="float32", order="F")
    )
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


def test_gene_fraction_regularization() -> None:
    """
    The regularization reaches the axes it is applied on.
    """
    graph = mg.metacells_gene_gene_graph(_test_daf(), x_gene="A", y_gene="B", gene_fraction_regularization=1e-3)
    assert graph.configuration.x_axis.log_regularization == 1e-3
    assert graph.configuration.y_axis.log_regularization == 1e-3
