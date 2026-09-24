"""
Test the data sources.

The Julia side tests what each source computes; what is worth testing here is that each kind of wrapper reaches it:
the sinks, the repository, the frame and the arrays all cross to Julia, and the results come back as Python values.
"""

from typing import Any
from typing import List

import dafpy as dp
import numpy as np
import pandas as pd
import somegraphspy as sg

import metacellsgraphspy as mg


def _list(values: Any) -> List[Any]:
    """
    A graph field as a list, once it is known to be filled.
    """
    assert values is not None
    return list(values)


def _test_daf() -> dp.DafWriter:
    """
    A repository with types, blocks and counts.
    """
    daf = dp.memory_daf(name="test!")
    daf.add_axis("gene", ["A", "B", "C"])
    daf.add_axis("metacell", ["M1", "M2", "M3"])
    daf.add_axis("block", ["B1", "B2"])
    daf.add_axis("type", ["T1", "T2"])
    daf.set_vector("metacell", "block", np.array(["B1", "B1", "B2"], dtype=str))
    daf.set_vector("metacell", "type", np.array(["T1", "T2", ""], dtype=str))
    daf.set_vector("metacell", "n_cells", np.array([10, 20, 30], dtype="uint32"))
    daf.set_vector("metacell", "total_UMIs", np.array([1000, 2000, 3000], dtype="uint64"))
    daf.set_vector("block", "total_UMIs", np.array([3000, 6000], dtype="uint64"))
    daf.set_vector("type", "color", np.array(["red", "blue"], dtype=str))
    daf.set_vector("gene", "marker_rank", np.array([2, 1, np.iinfo(np.uint32).max], dtype="uint32"))
    daf.set_vector("gene", "is_skeleton", np.array([False, True, False]))
    # Column-major, which is what ``Daf`` stores: the layout is not a detail it papers over.
    daf.set_matrix(
        "gene",
        "metacell",
        "linear_fraction",
        np.array([[0.1, 0.2, 0.3], [0.3, 0.2, 0.1], [0.5, 0.5, 0.5]], dtype="float32", order="F"),
    )
    return daf


def test_constants() -> None:
    """
    The constants are read from Julia.
    """
    assert mg.GENE_FRACTION_REGULARIZATION_FOR_GRAPHS == 1e-5
    assert mg.MAX_FOLD_FOR_GRAPHS == 3.0
    assert mg.EMPTY_TYPE_COLOR == "magenta"


def test_fill_a_points_graph() -> None:
    """
    Filling the views of a graph from a repository, the way the Julia graph functions do.
    """
    daf = _test_daf()
    graph = sg.points_graph()
    mg.fill_gene_expression(graph.x_axis_vector_fields(), daf, gene="A")
    mg.fill_gene_expression(graph.y_axis_vector_fields(), daf, gene="B")
    mg.fill_type(graph.points_colors_vector_fields(), daf)
    graph.validate()

    assert _list(graph.data.x.vector) == [np.float32(0.1), np.float32(0.2), np.float32(0.3)]
    assert graph.configuration.x_axis.title == "A fraction"
    assert graph.configuration.x_axis.scale.log_base == sg.LogBase.Log2Base
    assert _list(graph.data.points.entities.names) == ["M1", "M2", "M3"]
    assert _list(graph.data.points.colors.vector) == ["T1", "T2", ""]
    assert graph.configuration.points.colors.palette == {"T1": "red", "T2": "blue", "": "magenta"}


def test_fill_several_sinks() -> None:
    """
    One source feeds several places in the graph at once, and a mask hides entities.
    """
    daf = _test_daf()
    graph = sg.points_graph()
    mg.fill_n_cells([graph.x_axis_vector_fields(), graph.points_colors_vector_fields()], daf)
    mg.put_vector_mask_data(graph.points_entities(), [True, False, True])
    assert _list(graph.data.x.vector) == [10, 20, 30]
    assert _list(graph.data.points.colors.vector) == [10, 20, 30]
    assert graph.configuration.x_axis.scale.log_base == sg.LogBase.Log2Base
    assert graph.configuration.points.colors.scale.log_base == sg.LogBase.Log2Base
    assert _list(graph.data.points.entities.mask) == [True, False, True]


def test_get_vectors() -> None:
    """
    The vector getters return ``numpy`` arrays, whatever the data type.
    """
    daf = _test_daf()
    assert list(mg.get_axis_entries_vector(daf, axis="metacell", entries=[3, 1])) == ["M3", "M1"]
    assert list(mg.get_axis_vector(daf, axis="metacell", query_suffix=": block")) == ["B1", "B1", "B2"]
    assert list(mg.get_block_vector(daf)) == ["B1", "B1", "B2"]
    assert list(mg.get_type_vector(daf, entries=["M1", "M3"])) == ["T1", ""]
    assert mg.get_type_colors(daf) == {"T1": "red", "T2": "blue", "": "magenta"}
    assert list(mg.get_n_cells_vector(daf)) == [10, 20, 30]
    assert list(mg.get_total_UMIs_vector(daf, via=["block"])) == [3000, 3000, 6000]
    assert list(mg.get_top_marker_gene_indices(daf, markers_count=1)) == [2]
    assert list(mg.get_skeleton_gene_indices(daf)) == [2]
    assert list(mg.get_boolean_annotation_vector(daf, axis="gene", property="is_skeleton")) == [
        "false",
        "true",
        "false",
    ]
    assert list(mg.get_vector_query(daf, query="@ metacell : n_cells", indices=[2])) == [20]


def test_get_matrices() -> None:
    """
    The matrix getters return ``numpy`` matrices, with the rows and columns as Julia has them.
    """
    daf = _test_daf()
    matrix = mg.get_genes_expression_matrix(daf, genes=["A", "C"], entries=["M1"])
    assert matrix.shape == (2, 1)
    assert list(matrix[:, 0]) == [np.float32(0.1), np.float32(0.5)]

    matrix = mg.get_axes_matrix(daf, rows_axis="gene", columns_axis="metacell", query_suffix=":: linear_fraction")
    assert matrix.shape == (3, 3)

    matrix = mg.get_matrix_query(daf, query="@ gene @ metacell :: linear_fraction", row_indices=[3])
    assert matrix.shape == (1, 3)


def test_fill_a_heatmap_graph() -> None:
    """
    A matrix source names both axes, and an annotation of one axis shares its entities.
    """
    daf = _test_daf()
    graph = sg.heatmap_graph()
    mg.fill_genes_expression_matrix(graph.entries_matrix_fields(), daf, genes=[1, 2])
    mg.fill_type(graph.columns_annotations_colors_vector_fields(graph.add_columns_annotation()), daf)
    graph.validate()
    assert _list(graph.data.rows.entities.names) == ["A", "B"]
    assert _list(graph.data.columns.entities.names) == ["M1", "M2", "M3"]
    matrix = graph.data.entries.matrix
    assert matrix is not None
    assert matrix.shape == (2, 3)
    assert _list(graph.data.columns.annotations[0].values.vector) == ["T1", "T2", ""]


def test_frame_sources() -> None:
    """
    A ``pandas`` frame is a source too, the way a gene report computed by ``metacellspy`` is.
    """
    frame = pd.DataFrame({"gene": ["A", "B"], "score": [0.5, 1.5], "is_marker": [True, False]})
    graph = sg.bars_graph()
    mg.fill_column_vector_data(graph.values_axis_vector_fields(), frame, column="score", title="Score")
    mg.fill_column_names_data(graph.bars_entities(), frame, column="gene")
    mg.fill_column_boolean_annotation(
        graph.annotations_colors_vector_fields(graph.add_annotation()), frame, column="is_marker"
    )
    graph.validate()
    assert _list(graph.data.values.vector) == [0.5, 1.5]
    assert _list(graph.data.bars.hovers) == ["Score: 0.5<br>is_marker: true", "Score: 1.5<br>is_marker: false"]
    assert _list(graph.data.bars.names) == ["A", "B"]
    assert _list(graph.data.annotations[0].values.vector) == ["true", "false"]
    assert list(mg.get_column_vector(frame, column="score")) == [0.5, 1.5]


def test_put_data_and_configuration() -> None:
    """
    Data and configuration already in hand are written by the ``put_`` functions.
    """
    graph = sg.points_graph()
    mg.put_vector_data(graph.x_axis_vector_fields(), np.array([1.0, 2.0]), title="X")
    mg.put_vector_names_data(graph.points_entities(), ["a", "b"])
    mg.put_umap_data(graph.y_axis_vector_fields(), np.array([3.0, 4.0]))
    mg.put_umap_configuration(graph.y_axis_vector_fields())
    mg.put_count_configuration(graph.x_axis_vector_fields())
    mg.put_type_configuration(graph.points_colors_vector_fields(), {"a": "red"}, title="kind")
    assert _list(graph.data.x.vector) == [1.0, 2.0]
    assert _list(graph.data.points.entities.names) == ["a", "b"]
    assert _list(graph.data.y.vector) == [3.0, 4.0]
    assert not graph.configuration.y_axis.show_ticks
    assert graph.configuration.x_axis.scale.log_base == sg.LogBase.Log2Base
    assert graph.configuration.points.colors.palette == {"a": "red"}
    assert graph.configuration.points.colors.title == "kind"
