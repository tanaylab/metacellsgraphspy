"""
Test the bar graphs.

The Julia side tests what the graphs contain; what is worth testing here is that the wrapper reaches it and hands back
something the notebook can display - which is a ``somegraphspy`` graph rather than a raw Julia object.
"""

# The base repository is built the same way as in ``test_scatter_graphs.py``, since both graphs are read against
# the same shape of base. The same disable is in that file, since the check needs both sides to report it.
# pylint: disable=duplicate-code

import dafpy as dp
import numpy as np
from somegraphspy import SeriesBarsGraph

import metacellsgraphspy as mg

CORRELATION = "correlation_between_base_neighborhood_cells_and_punctuated_metacells"


def _test_dafs() -> tuple[dp.DafWriter, dp.DafWriter]:
    """
    A base repository and a repository scored against it, with a known count of improved and degraded blocks per gene.
    """
    base_daf = dp.memory_daf(name="base!")
    base_daf.add_axis("gene", ["A", "B"])
    base_daf.add_axis("block", ["B1", "B2"])
    base_daf.add_axis("base_block", ["B1", "B2"])
    # Column-major, which is what ``Daf`` stores: the layout is not a detail it papers over.
    base_daf.set_matrix(
        "gene", "base_block", CORRELATION, np.array([[0.5, 0.5], [0.5, 0.5]], dtype="float32", order="F")
    )

    daf = dp.memory_daf(name="test!")
    daf.add_axis("gene", ["A", "B"])
    daf.add_axis("base_block", ["B1", "B2"])
    # The graphs draw the gene report, which is over the marker genes and says which of them are what.
    daf.set_vector("gene", "is_marker", np.array([True, True]))
    daf.set_vector("gene", "marker_rank", np.array([1, 2], dtype="uint32"))
    daf.set_vector("gene", "is_lateral", np.array([False, True]))
    daf.set_vector("gene", "is_transcription_factor", np.array([True, False]))
    daf.set_vector("gene", "is_regulator", np.array([True, False]))
    daf.set_vector("gene", "is_skeleton", np.array([False, True]))
    daf.set_matrix(
        "gene", "gene", "correlation_between_markers", np.array([[1.0, 0.5], [0.5, 1.0]], dtype="float32", order="F")
    )
    # ``A`` improves in both blocks and ``B`` degrades in both.
    daf.set_matrix("gene", "base_block", CORRELATION, np.array([[0.6, 0.6], [0.4, 0.4]], dtype="float32", order="F"))
    # Nothing shares a module with anything, which is enough for the hovers to have something to say.
    for side in ("improved", "degraded"):
        daf.set_vector(
            "gene",
            f"mean_no_module_fraction_in_base_neighborhood_cells_at_{side}_base_blocks",
            np.array([1.0, 1.0], dtype="float32"),
        )
        daf.set_matrix(
            "gene",
            "gene",
            f"mean_shared_module_fraction_in_base_neighborhood_cells_at_{side}_base_blocks",
            np.array([[0.0, 0.0], [0.0, 0.0]], dtype="float32", order="F"),
        )

    return (daf, base_daf)


def test_improved_genes_graph() -> None:
    """
    The improved genes graph is a mirrored series bars graph, best gene last.
    """
    daf, base_daf = _test_dafs()
    graph = mg.improved_genes_graph(daf=daf, base_daf=base_daf)
    assert isinstance(graph, SeriesBarsGraph)
    bars_names = graph.data.bars.names
    assert bars_names is not None
    assert list(bars_names) == ["B", "A"]
    assert graph.configuration.mirrored
    assert graph.figure is not None


def test_degraded_genes_graph() -> None:
    """
    The degraded genes graph picks its genes by the other side.
    """
    daf, base_daf = _test_dafs()
    graph = mg.degraded_genes_graph(daf=daf, base_daf=base_daf)
    assert isinstance(graph, SeriesBarsGraph)
    bars_names = graph.data.bars.names
    assert bars_names is not None
    assert list(bars_names) == ["A", "B"]


def test_genes_count() -> None:
    """
    The count reaches the choice of which genes are shown.
    """
    daf, base_daf = _test_dafs()
    graph = mg.improved_genes_graph(daf=daf, base_daf=base_daf, genes_count=1)
    bars_names = graph.data.bars.names
    assert bars_names is not None
    assert list(bars_names) == ["A"]
