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
    assert mg.EMPTY_CATEGORY_COLOR == "grey"


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
    assert graph.configuration.points.colors.palette == {"T1": "red", "T2": "blue", "": "grey"}


def test_fill_several_sinks() -> None:
    """
    One source feeds several places in the graph at once, and a mask hides entities.
    """
    daf = _test_daf()
    graph = sg.points_graph()
    mg.fill_n_cells([graph.x_axis_vector_fields(), graph.points_colors_vector_fields()], daf)
    sg.put_vector_mask_data(graph.points_entities(), [True, False, True])
    sg.put_vector_order_data(graph.points_entities(), [3, 1, 2])
    assert _list(graph.data.x.vector) == [10, 20, 30]
    assert _list(graph.data.points.colors.vector) == [10, 20, 30]
    assert graph.configuration.x_axis.scale.log_base == sg.LogBase.Log10Base
    assert graph.configuration.points.colors.scale.log_base == sg.LogBase.Log10Base
    assert _list(graph.data.points.entities.mask) == [True, False, True]
    assert _list(graph.data.points.entities.order) == [3, 1, 2]


def test_get_vectors() -> None:
    """
    The vector getters return ``numpy`` arrays, whatever the data type.
    """
    daf = _test_daf()
    assert list(mg.get_axis_entries_vector(daf, axis="metacell", entries=[3, 1])) == ["M3", "M1"]
    assert list(mg.get_axis_vector(daf, axis="metacell", query_suffix=": block")) == ["B1", "B1", "B2"]
    assert list(mg.get_block_vector(daf)) == ["B1", "B1", "B2"]
    assert list(mg.get_block_colors(daf)) == ["B1", "B2", ""]
    assert mg.get_block_colors(daf)[""] == "grey"
    assert list(mg.get_type_vector(daf, entries=["M1", "M3"])) == ["T1", ""]
    assert mg.get_type_colors(daf) == {"T1": "red", "T2": "blue", "": "grey"}
    assert list(mg.get_n_cells_vector(daf)) == [10, 20, 30]
    assert list(mg.get_total_UMIs_vector(daf, via=["block"])) == [3000, 3000, 6000]
    assert list(mg.get_top_marker_gene_indices(daf, markers_count=1)) == [2]
    assert list(mg.get_skeleton_gene_indices(daf)) == [2]
    assert list(mg.get_boolean_vector(daf, axis="gene", property="is_skeleton")) == [False, True, False]
    assert list(mg.get_is_skeleton_vector(daf)) == [False, True, False]
    assert list(mg.get_axis_entries_vector(daf, axis="gene", entries=mg.get_is_skeleton_vector(daf))) == ["B"]
    assert list(mg.get_size_vector(daf, axis="metacell", property="n_cells")) == [10, 20, 30]
    assert list(mg.get_category_vector(daf, axis="metacell", property="type")) == ["T1", "T2", ""]
    assert mg.get_category_colors(daf, axis="type") == mg.get_type_colors(daf)
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
    mg.fill_block(graph.columns_annotations_colors_vector_fields(graph.add_columns_annotation()), daf)
    graph.validate()
    assert _list(graph.data.rows.entities.names) == ["A", "B"]
    assert _list(graph.data.columns.entities.names) == ["M1", "M2", "M3"]
    matrix = graph.data.entries.matrix
    assert matrix is not None
    assert matrix.shape == (2, 3)
    assert _list(graph.data.columns.annotations[0].values.vector) == ["T1", "T2", ""]
    assert _list(graph.data.columns.annotations[1].values.vector) == ["B1", "B1", "B2"]


def test_fill_named_entries() -> None:
    """
    A fill which isn't given its entries fills the entries its sinks are already named after.
    """
    daf = _test_daf()
    graph = sg.heatmap_graph()
    mg.fill_genes_expression_matrix(graph.entries_matrix_fields(), daf, genes=["A"], entries=["M3", "M1"])
    mg.fill_block(graph.columns_annotations_colors_vector_fields(graph.add_columns_annotation()), daf)
    graph.validate()
    assert _list(graph.data.columns.annotations[0].values.vector) == ["B2", "B1"]
    assert _list(graph.data.columns.entities.names) == ["M3", "M1"]


def test_bucket_fills() -> None:
    """
    The generic fill and put of each display bucket, and a per-property fill, reach Julia.
    """
    daf = _test_daf()
    graph = sg.points_graph()
    mg.fill_size(graph.x_axis_vector_fields(), daf, axis="metacell", property="n_cells", log_base=sg.LogBase.Log2Base)
    mg.fill_category(graph.points_colors_vector_fields(), daf, axis="metacell", property="type")
    assert _list(graph.data.x.vector) == [10, 20, 30]
    assert graph.configuration.x_axis.scale.log_base == sg.LogBase.Log2Base
    assert _list(graph.data.points.colors.vector) == ["T1", "T2", ""]

    graph = sg.points_graph()
    names = ["M1", "M2", "M3"]
    mg.put_size_values(graph.x_axis_vector_fields(), np.array([1.0, 2.0, 3.0]), names)
    mg.put_category_values(
        graph.points_colors_vector_fields(), np.array(["B1", "B1", "B2"]), names, mg.get_block_colors(daf)
    )
    assert _list(graph.data.x.vector) == [1.0, 2.0, 3.0]
    assert graph.configuration.points.colors.palette == mg.get_block_colors(daf)

    heatmap = sg.heatmap_graph()
    mg.fill_genes_expression_matrix(heatmap.entries_matrix_fields(), daf, genes=["A", "B"])
    mg.fill_is_skeleton(heatmap.rows_annotations_colors_vector_fields(heatmap.add_rows_annotation()), daf)
    mg.fill_boolean(
        heatmap.rows_annotations_colors_vector_fields(heatmap.add_rows_annotation()),
        daf,
        axis="gene",
        property="is_skeleton",
    )
    mg.put_boolean_values(
        heatmap.rows_annotations_colors_vector_fields(heatmap.add_rows_annotation()),
        np.array([True, False]),
        ["A", "B"],
    )
    heatmap.validate()
    assert _list(heatmap.data.rows.annotations[0].values.vector) == ["false", "true"]
    assert _list(heatmap.data.rows.annotations[1].values.vector) == ["false", "true"]
    assert _list(heatmap.data.rows.annotations[2].values.vector) == ["true", "false"]

    daf.set_vector(
        "gene", "correlation_between_cells_and_punctuated_metacells", np.array([0.5, -0.1, 0.2], dtype="float32")
    )
    graph = sg.points_graph()
    mg.fill_correlation_between_cells_and_punctuated_metacells(graph.points_colors_vector_fields(), daf)
    assert graph.configuration.points.colors.scale.minimum == np.float32(-0.1)
    mg.fill_correlation(
        graph.points_colors_vector_fields(),
        daf,
        axis="gene",
        property="correlation_between_cells_and_punctuated_metacells",
    )
    mg.put_correlation_values(graph.points_colors_vector_fields(), np.array([0.1, 0.2, 0.3]), ["A", "B", "C"])
    assert graph.configuration.points.colors.scale.minimum == 0
    assert list(mg.get_correlation_between_cells_and_punctuated_metacells_vector(daf)) == list(
        mg.get_correlation_vector(daf, axis="gene", property="correlation_between_cells_and_punctuated_metacells")
    )


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


def test_hover_format() -> None:
    """
    A generic fill passes its ``hover_format`` to the hover line, leaving the values as they are.
    """
    frame = pd.DataFrame({"share": [0.25, 1.0]})
    graph = sg.bars_graph()
    mg.fill_column_vector_data(
        graph.values_axis_vector_fields(),
        frame,
        column="share",
        title="Share",
        hover_format=sg.percent_format(digits=1),
    )
    assert _list(graph.data.values.vector) == [0.25, 1.0]
    assert _list(graph.data.bars.hovers) == ["Share: 25.0%", "Share: 100.0%"]

    daf = _test_daf()
    points = sg.points_graph()
    mg.fill_axis_vector_data(
        points.x_axis_vector_fields(),
        daf,
        axis="metacell",
        query_suffix=": n_cells",
        title="Cells",
        hover_format=sg.scientific_format(digits=1),
    )
    assert _list(points.data.points.entities.hovers) == ["Cells: 1.0e+01", "Cells: 2.0e+01", "Cells: 3.0e+01"]

    heatmap = sg.heatmap_graph()
    mg.fill_axes_matrix_data(
        heatmap.entries_matrix_fields(),
        daf,
        rows_axis="gene",
        columns_axis="metacell",
        query_suffix=":: linear_fraction",
        hover_format=sg.fixed_format(digits=2),
    )
    cells_hovers = heatmap.data.cells.hovers
    assert cells_hovers is not None
    assert cells_hovers[0, 0] == "0.10"


def test_bucket_hover_formats() -> None:
    """
    The data sources of a bucket show the hovers by the bucket's hover format, unless given another one.
    """
    daf = _test_daf()
    points = sg.points_graph()
    mg.fill_n_cells(points.x_axis_vector_fields(), daf)
    assert _list(points.data.points.entities.hovers) == ["cells: 10", "cells: 20", "cells: 30"]

    points = sg.points_graph()
    mg.put_correlation_values(points.x_axis_vector_fields(), np.array([0.1234, -0.5]), ["A", "B"], title="Correlation")
    assert _list(points.data.points.entities.hovers) == ["Correlation: 0.123", "Correlation: -0.500"]

    points = sg.points_graph()
    mg.put_size_values(
        points.x_axis_vector_fields(),
        np.array([12.34, 5.0]),
        ["A", "B"],
        title="Mean",
        hover_format=mg.MEAN_SIZE_HOVER_FORMAT_FOR_GRAPHS,
    )
    assert _list(points.data.points.entities.hovers) == ["Mean: 12.3", "Mean: 5.0"]


def test_put_data_and_configuration() -> None:
    """
    Data and configuration already in hand are written by the ``put_`` functions.
    """
    graph = sg.points_graph()
    sg.put_vector_data(graph.x_axis_vector_fields(), np.array([1.0, 2.0]), title="X")
    sg.put_vector_names_data(graph.points_entities(), ["a", "b"])
    mg.put_umap_data(graph.y_axis_vector_fields(), np.array([3.0, 4.0]))
    mg.put_coordinates_configuration(graph.y_axis_vector_fields())
    mg.put_size_configuration(graph.x_axis_vector_fields())
    mg.put_category_configuration(graph.points_colors_vector_fields(), {"a": "red"}, title="kind", show_legend=False)
    assert _list(graph.data.x.vector) == [1.0, 2.0]
    assert _list(graph.data.points.entities.names) == ["a", "b"]
    assert _list(graph.data.y.vector) == [3.0, 4.0]
    assert not graph.configuration.y_axis.show_ticks
    assert graph.configuration.x_axis.scale.log_base == sg.LogBase.Log10Base
    assert graph.configuration.x_axis.scale.log_regularization == 1
    assert graph.configuration.points.colors.palette == {"a": "red"}
    assert graph.configuration.points.colors.title == "kind"
    assert not graph.configuration.points.colors.show_legend

    graph = sg.points_graph()
    mg.put_correlation_configuration(graph.points_colors_vector_fields(), np.array([0.2, 0.6]))
    assert graph.configuration.points.colors.scale.minimum == 0
    assert graph.configuration.points.colors.scale.maximum == 0.6


def test_entry_masks() -> None:
    """
    The entry masks cover the whole axis, and combine to pick the entries of another data source.
    """
    daf = _test_daf()
    daf.set_matrix("block", "block", "is_in_neighborhood", np.array([[True, True], [False, True]], order="F"))
    daf.set_matrix(
        "metacell", "block", "is_in_environment", np.array([[True, False], [True, True], [False, True]], order="F")
    )
    is_t2 = mg.get_is_in_categories_vector(daf, axis="metacell", property="type", categories=["T2"])
    assert list(is_t2) == [False, True, False]
    assert list(mg.get_is_in_neighborhood_vector(daf, block="B1")) == [True, True, False]
    assert list(mg.get_is_in_neighborhood_vector(daf, axis="block", block="B2")) == [True, True]
    is_in_environment = mg.get_is_in_environment_vector(daf, block="B2")
    assert list(is_in_environment) == [False, True, True]
    assert list(mg.get_axis_entries_vector(daf, axis="metacell", entries=is_t2 & is_in_environment)) == ["M2"]
