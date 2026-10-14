"""
Test the ``AnnData`` graphs.
"""

from pathlib import Path
from typing import Any
from typing import List

import numpy as np
import pandas as pd
import somegraphspy as sg
from anndata import AnnData  # type: ignore
from somegraphspy import HeatmapGraph
from somegraphspy import PointsGraph

import metacellsgraphspy as mg


def _list(values: Any) -> List[Any]:
    """
    A graph field as a list, once it is known to be filled.
    """
    assert values is not None
    return list(values)


def _test_adata() -> AnnData:
    """
    Metacells with fractions of two genes, and a type for all but one of them.
    """
    return AnnData(
        X=np.array([[0.1, 0.3], [0.2, 0.2], [0.3, 0.1]], dtype="float32"),
        obs=pd.DataFrame({"cell_type": pd.Categorical(["T1", None, "T2"])}, index=["M1", "M2", "M3"]),
        var=pd.DataFrame(index=["A", "B"]),
    )


def test_ad_gene_gene_graph() -> None:
    """
    The gene-gene graph is a points graph with a point per metacell, the same as the ``Daf`` graph.
    """
    graph = mg.ad_gene_gene_graph(_test_adata(), x_gene="A", y_gene="B")
    assert isinstance(graph, PointsGraph)
    assert _list(graph.data.x.vector) == [np.float32(0.1), np.float32(0.2), np.float32(0.3)]
    assert _list(graph.data.y.vector) == [np.float32(0.3), np.float32(0.2), np.float32(0.1)]
    assert _list(graph.data.points.entities.names) == ["M1", "M2", "M3"]
    assert graph.configuration.y_axis.title == "B fraction"
    assert graph.figure is not None


def test_ad_gene_gene_graph_of_entries() -> None:
    """
    The entries pick the points.
    """
    graph = mg.ad_gene_gene_graph(_test_adata(), x_gene="A", y_gene="B", entries=["M3", "M1"])
    assert _list(graph.data.x.vector) == [np.float32(0.3), np.float32(0.1)]
    assert _list(graph.data.points.entities.names) == ["M3", "M1"]


def test_ad_gene_gene_graph_by_type(tmp_path: Path) -> None:
    """
    Coloring the points by type gives a graph which is valid and renders.
    """
    type_colors_csv = tmp_path / "type_colors.csv"
    type_colors_csv.write_text("type,color\nT1,#ff0000\nT2,#00ff00\n")
    adata = _test_adata()
    graph = mg.ad_gene_gene_graph(adata, x_gene="A", y_gene="B")
    mg.ad_fill_type(
        graph.points_colors_vector_fields(), adata, type_property="cell_type", type_colors_csv=str(type_colors_csv)
    )
    graph.validate()
    assert graph.figure is not None


def test_ad_genes_heatmap_graph() -> None:
    """
    The genes heatmap is a heatmap with a row per gene (by default, the strongest markers) and a column per metacell,
    laid out as the ``Daf`` graph is.
    """
    adata = _test_adata()
    adata.var["marker_gene"] = [True, True]
    graph = mg.ad_genes_heatmap_graph(adata)
    assert isinstance(graph, HeatmapGraph)
    assert _list(graph.data.rows.entities.names) == ["A", "B"]
    assert _list(graph.data.columns.entities.names) == ["M1", "M2", "M3"]
    assert graph.configuration.rows.title == "Genes"
    assert graph.configuration.columns.title == "Metacells"
    assert not graph.configuration.columns.show_ticks
    assert graph.figure is not None


def test_ad_markers_heatmap_graph() -> None:
    """
    A markers heatmap picks its genes by the marker rank, computed when the ``AnnData`` has none.
    """
    adata = _test_adata()
    adata.var["marker_gene"] = [True, True]
    graph = mg.ad_genes_heatmap_graph(
        adata, genes=mg.ad_get_top_marker_gene_indices(adata, markers_count=1), rows_axis_title="Markers"
    )
    assert _list(graph.data.rows.entities.names) == ["B"]
    assert graph.configuration.rows.title == "Markers"
    assert graph.configuration.rows.order_source is None
    assert graph.configuration.columns.order_source == sg.OrderSource.OptimalTreeReorder


def test_ad_default_annotations(tmp_path: Path) -> None:
    """
    The default annotations of each side are added when the ``AnnData`` has their columns.
    """
    type_colors_csv = tmp_path / "type_colors.csv"
    type_colors_csv.write_text("type,color\nT1,#ff0000\nT2,#00ff00\n")
    adata = _test_adata()
    graph = mg.ad_genes_heatmap_graph(adata, genes=["A", "B"])
    mg.ad_add_genes_heatmap_default_rows_annotations(graph, adata)
    mg.ad_add_genes_heatmap_default_columns_annotations(
        graph, adata, type_property="missing", type_colors_csv=str(type_colors_csv)
    )
    assert not graph.data.rows.annotations
    assert not graph.data.columns.annotations

    adata.var["lateral_gene"] = [False, True]
    mg.ad_add_genes_heatmap_default_rows_annotations(graph, adata)
    mg.ad_add_genes_heatmap_default_columns_annotations(
        graph, adata, type_property="cell_type", type_colors_csv=str(type_colors_csv)
    )
    assert _list(graph.data.rows.annotations[0].values.vector) == ["false", "true"]
    assert _list(graph.data.columns.annotations[0].values.vector) == ["T1", "", "T2"]
