"""
Test the bar graphs.

The Julia side tests what the graphs contain; what is worth testing here is that the wrapper reaches it and hands back
something the notebook can display - which is a ``somegraphspy`` graph rather than a raw Julia object.
"""

import dafpy as dp
import numpy as np
from somegraphspy import SeriesBarsGraph

import metacellsgraphspy as mg

CORRELATION = "correlation_between_base_neighborhood_cells_and_punctuated_metacells"


def _test_dafs() -> tuple[dp.DafWriter, dp.DafWriter]:
    """
    A base repository and a repository scored against it, with a known count of improved and declined blocks per gene.
    """
    base_daf = dp.memory_daf(name="base!")
    base_daf.add_axis("gene", ["A", "B"])
    base_daf.add_axis("block", ["B1", "B2"])
    base_daf.add_axis("base_block", ["B1", "B2"])
    base_daf.set_vector("gene", "is_lateral", np.array([False, True]))
    base_daf.set_vector("gene", "is_regulator", np.array([True, False]))
    # Column-major, which is what ``Daf`` stores: the layout is not a detail it papers over.
    base_daf.set_matrix(
        "gene", "base_block", CORRELATION, np.array([[0.5, 0.5], [0.5, 0.5]], dtype="float32", order="F")
    )

    daf = dp.memory_daf(name="test!")
    daf.add_axis("gene", ["A", "B"])
    daf.add_axis("base_block", ["B1", "B2"])
    # ``A`` improves in both blocks and ``B`` declines in both.
    daf.set_matrix("gene", "base_block", CORRELATION, np.array([[0.6, 0.6], [0.4, 0.4]], dtype="float32", order="F"))

    return (daf, base_daf)


def test_improved_genes_graph() -> None:
    """
    The improved genes graph is a mirrored series bars graph, best gene last.
    """
    daf, base_daf = _test_dafs()
    graph = mg.improved_genes_graph(daf=daf, base_daf=base_daf)
    assert isinstance(graph, SeriesBarsGraph)
    bars_names = graph.data.bars_names
    assert bars_names is not None
    assert list(bars_names) == ["B", "A"]
    assert graph.configuration.mirrored
    assert graph.figure is not None


def test_declined_genes_graph() -> None:
    """
    The declined genes graph picks its genes by the other side.
    """
    daf, base_daf = _test_dafs()
    graph = mg.declined_genes_graph(daf=daf, base_daf=base_daf)
    assert isinstance(graph, SeriesBarsGraph)
    bars_names = graph.data.bars_names
    assert bars_names is not None
    assert list(bars_names) == ["A", "B"]


def test_genes_count() -> None:
    """
    The count reaches the choice of which genes are shown.
    """
    daf, base_daf = _test_dafs()
    graph = mg.improved_genes_graph(daf=daf, base_daf=base_daf, genes_count=1)
    bars_names = graph.data.bars_names
    assert bars_names is not None
    assert list(bars_names) == ["A"]
