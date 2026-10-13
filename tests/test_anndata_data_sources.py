"""
Test the ``AnnData`` data sources.
"""

from pathlib import Path
from typing import Any
from typing import List
from typing import Mapping

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp  # type: ignore
import somegraphspy as sg
from anndata import AnnData  # type: ignore

import metacellsgraphspy as mg


def _list(values: Any) -> List[Any]:
    """
    A graph field as a list, once it is known to be filled.
    """
    assert values is not None
    return list(values)


def _test_adata() -> AnnData:
    """
    Metacells with fractions of two genes, a type for all but one of them, a UMAP and total UMIs, and the gene flags.
    """
    adata = AnnData(
        X=np.array([[0.1, 0.3], [0.2, 0.2], [0.3, 0.1]], dtype="float32"),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(["T1", None, "T2"]),
                "x": [0.0, 1.0, 2.0],
                "y": [5.0, 4.0, 3.0],
                "total_umis": [100.0, 200.0, 300.0],
            },
            index=["M1", "M2", "M3"],
        ),
        var=pd.DataFrame(
            {"lateral_gene": [True, False], "excluded_gene": [False, False], "marker_gene": [True, True]},
            index=["A", "B"],
        ),
    )
    return adata


def _type_colors_csv(tmp_path: Path) -> str:
    """
    A type colors file listing a type nothing has, in an order which isn't sorted.
    """
    path = tmp_path / "type_colors.csv"
    path.write_text('"cell_type","color"\n"T2","#00ff00"\n"T1","#ff0000"\n"T3","#0000ff"\n')
    return str(path)


def test_ad_get_gene_expression_vector() -> None:
    """
    The vector is the gene's column of ``X``, for the entries.
    """
    adata = _test_adata()
    assert _list(mg.ad_get_gene_expression_vector(adata, gene="B")) == [
        np.float32(0.3),
        np.float32(0.2),
        np.float32(0.1),
    ]
    assert _list(mg.ad_get_gene_expression_vector(adata, gene="A", entries=["M3", "M1"])) == [
        np.float32(0.3),
        np.float32(0.1),
    ]
    assert _list(mg.ad_get_gene_expression_vector(adata, gene="A", entries=[2])) == [np.float32(0.2)]


def test_ad_get_gene_expression_vector_of_sparse() -> None:
    """
    A sparse ``X`` gives the same vector as a dense one.
    """
    adata = _test_adata()
    adata.X = sp.csr_matrix(adata.X)
    assert _list(mg.ad_get_gene_expression_vector(adata, gene="B")) == [
        np.float32(0.3),
        np.float32(0.2),
        np.float32(0.1),
    ]


def test_ad_unknown_entries() -> None:
    """
    An entry name which isn't in the ``obs`` is an error.
    """
    with pytest.raises(KeyError, match="M4"):
        mg.ad_get_gene_expression_vector(_test_adata(), gene="A", entries=["M1", "M4"])


def test_ad_masked_entries() -> None:
    """
    A Boolean mask picks the entries, and must cover the whole ``obs``.
    """
    adata = _test_adata()
    assert _list(mg.ad_get_gene_expression_vector(adata, gene="A", entries=np.array([True, False, True]))) == [
        np.float32(0.1),
        np.float32(0.3),
    ]
    with pytest.raises(ValueError, match="mask length: 2"):
        mg.ad_get_gene_expression_vector(adata, gene="A", entries=[True, False])


def test_ad_fill_gene_expression() -> None:
    """
    The gene expression fills an axis the same way as the ``Daf`` data source.
    """
    graph = sg.points_graph()
    mg.ad_fill_gene_expression(graph.x_axis_vector_fields(), _test_adata(), gene="A")
    assert _list(graph.data.x.vector) == [np.float32(0.1), np.float32(0.2), np.float32(0.3)]
    assert all(hover.startswith("A fraction: ") for hover in _list(graph.data.points.entities.hovers))
    assert _list(graph.data.points.entities.names) == ["M1", "M2", "M3"]
    assert graph.configuration.x_axis.title == "A fraction"
    assert graph.configuration.x_axis.scale.log_base == sg.LogBase.Log2Base
    assert graph.configuration.x_axis.scale.log_regularization == mg.GENE_FRACTION_REGULARIZATION_FOR_GRAPHS


def test_ad_get_type_vector() -> None:
    """
    A missing type is the empty string.
    """
    adata = _test_adata()
    assert _list(mg.ad_get_type_vector(adata, type_property="cell_type")) == ["T1", "", "T2"]
    assert _list(mg.ad_get_type_vector(adata, type_property="cell_type", entries=[3, 2])) == ["T2", ""]


def test_ad_get_type_colors(tmp_path: Path) -> None:
    """
    The palette follows the file, and adds a color for the empty type.
    """
    type_colors_csv = _type_colors_csv(tmp_path)
    assert list(mg.ad_get_type_colors(type_colors_csv).items()) == [
        ("T2", "#00ff00"),
        ("T1", "#ff0000"),
        ("T3", "#0000ff"),
        ("", mg.EMPTY_CATEGORY_COLOR),
    ]
    assert mg.ad_get_type_colors(type_colors_csv, empty_type_color="grey")[""] == "grey"


def test_ad_fill_type(tmp_path: Path) -> None:
    """
    The types fill the points colors the same way as the ``Daf`` data source.
    """
    graph = sg.points_graph()
    mg.ad_fill_type(
        graph.points_colors_vector_fields(),
        _test_adata(),
        type_property="cell_type",
        type_colors_csv=_type_colors_csv(tmp_path),
    )
    assert _list(graph.data.points.colors.vector) == ["T1", "", "T2"]
    assert _list(graph.data.points.entities.hovers) == ["type: T1", "type: ", "type: T2"]
    assert _list(graph.data.points.entities.names) == ["M1", "M2", "M3"]
    palette = graph.configuration.points.colors.palette
    assert isinstance(palette, Mapping)
    assert dict(palette) == {
        "T1": "#ff0000",
        "T2": "#00ff00",
        "T3": "#0000ff",
        "": mg.EMPTY_CATEGORY_COLOR,
    }
    assert graph.configuration.points.colors.show_legend
    assert graph.configuration.points.colors.title == "type"


def test_ad_fill_named_entries(tmp_path: Path) -> None:
    """
    A fill which isn't given its entries fills the entries its sinks are already named after.
    """
    graph = sg.points_graph()
    adata = _test_adata()
    mg.ad_fill_gene_expression(graph.x_axis_vector_fields(), adata, gene="A", entries=["M3", "M1"])
    mg.ad_fill_type(
        graph.points_colors_vector_fields(),
        adata,
        type_property="cell_type",
        type_colors_csv=_type_colors_csv(tmp_path),
    )
    assert _list(graph.data.points.colors.vector) == ["T2", "T1"]
    assert _list(graph.data.points.entities.names) == ["M3", "M1"]


def test_ad_umap() -> None:
    """
    The UMAP coordinates are the ``x`` and ``y`` of the ``obs``, shown without ticks.
    """
    adata = _test_adata()
    assert _list(mg.ad_get_umap_vector(adata, coordinate="y", entries=[3])) == [3.0]
    graph = sg.points_graph()
    mg.ad_fill_umap(graph.x_axis_vector_fields(), adata, coordinate="x")
    mg.ad_fill_umap(graph.y_axis_vector_fields(), adata, coordinate="y")
    assert _list(graph.data.x.vector) == [0.0, 1.0, 2.0]
    assert graph.configuration.x_axis.title == "UMAP X"
    assert not graph.configuration.x_axis.show_ticks
    assert not graph.configuration.y_axis.show_grid
    assert _list(graph.data.points.entities.names) == ["M1", "M2", "M3"]


def test_ad_total_UMIs() -> None:  # pylint: disable=invalid-name
    """
    The total UMIs are the ``total_umis`` of the ``obs``, shown as a size.
    """
    adata = _test_adata()
    assert _list(mg.ad_get_total_UMIs_vector(adata, entries=["M2"])) == [200.0]
    graph = sg.points_graph()
    mg.ad_fill_total_UMIs(graph.points_colors_vector_fields(), adata, log_base=sg.LogBase.Log2Base)
    assert _list(graph.data.points.colors.vector) == [100.0, 200.0, 300.0]
    assert graph.configuration.points.colors.title == "total UMIs"
    assert graph.configuration.points.colors.scale.log_base == sg.LogBase.Log2Base


def test_ad_gene_flags() -> None:
    """
    The gene flags are Boolean columns of the ``var``, whose entries are genes, and any of them serves as a mask.
    """
    adata = _test_adata()
    assert _list(mg.ad_get_is_lateral_vector(adata)) == [True, False]
    assert _list(mg.ad_get_is_excluded_vector(adata, entries=["B"])) == [False]
    assert _list(mg.ad_get_is_marker_vector(adata, entries=mg.ad_get_is_lateral_vector(adata))) == [True]

    heatmap = sg.heatmap_graph()
    mg.ad_fill_is_lateral(heatmap.rows_annotations_colors_vector_fields(heatmap.add_rows_annotation()), adata)
    mg.ad_fill_is_excluded(heatmap.rows_annotations_colors_vector_fields(heatmap.add_rows_annotation()), adata)
    mg.ad_fill_is_marker(heatmap.rows_annotations_colors_vector_fields(heatmap.add_rows_annotation()), adata)
    assert _list(heatmap.data.rows.entities.names) == ["A", "B"]
    assert _list(heatmap.data.rows.annotations[0].values.vector) == ["true", "false"]
    assert heatmap.data.rows.annotations[0].colors.title == "lateral"
    assert heatmap.data.rows.annotations[1].colors.title == "excluded"
    assert heatmap.data.rows.annotations[2].colors.title == "marker"


def test_ad_genes_expression_matrix() -> None:
    """
    The matrix is ``X`` as stored, with a row per gene and a column per metacell, from a dense or a sparse ``X``.
    """
    adata = _test_adata()
    matrix = mg.ad_get_genes_expression_matrix(adata, genes=["B"], entries=[1, 3])
    assert matrix.shape == (1, 2)
    assert list(matrix[0, :]) == [np.float32(0.3), np.float32(0.1)]

    adata.X = sp.csr_matrix(adata.X)
    heatmap = sg.heatmap_graph()
    mg.ad_fill_genes_expression_matrix(heatmap.entries_matrix_fields(), adata)
    entries_matrix = heatmap.data.entries.matrix
    assert entries_matrix is not None
    assert list(entries_matrix[:, 0]) == [np.float32(0.1), np.float32(0.3)]
    assert _list(heatmap.data.rows.entities.names) == ["A", "B"]
    assert _list(heatmap.data.columns.entities.names) == ["M1", "M2", "M3"]
    assert heatmap.configuration.entries.colors.title == "fraction"


def test_ad_genes_fold_matrix() -> None:
    """
    The fold is the log of the expression level, regularized as ``Metacells`` does, minus each gene's median.
    """
    adata = _test_adata()
    folds = mg.ad_get_genes_fold_matrix(adata)
    assert [round(float(fold), 3) for fold in folds[0, :]] == [-1.0, 0.0, 0.585]

    heatmap = sg.heatmap_graph()
    mg.ad_fill_genes_fold_matrix(heatmap.entries_matrix_fields(), adata, genes=["B"], entries=["M3", "M1"])
    matrix = heatmap.data.entries.matrix
    assert matrix is not None
    assert [round(float(fold), 3) for fold in matrix[0, :]] == [-0.792, 0.792]
    assert _list(heatmap.data.columns.entities.names) == ["M3", "M1"]
    assert heatmap.configuration.entries.colors.palette == "BuWtRd"
