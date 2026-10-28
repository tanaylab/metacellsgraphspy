"""
Make the calls of a typical session (for internal use).

The Julia package ``MetacellsGraphsPy`` runs this while it is precompiled, so the code these calls compile is cached, and
a session does not compile it again. It is not part of the API, so the package does not import it.
"""

import os
from tempfile import TemporaryDirectory

import dafpy as dp
import numpy as np
import pandas as pd
from anndata import AnnData  # type: ignore

from .anndata_data_sources import ad_fill_gene_expression
from .anndata_data_sources import ad_fill_is_excluded
from .anndata_data_sources import ad_fill_is_marker
from .anndata_data_sources import ad_fill_total_UMIs
from .anndata_data_sources import ad_fill_type
from .anndata_data_sources import ad_fill_umap
from .anndata_data_sources import ad_get_top_marker_gene_indices
from .anndata_graphs import ad_add_genes_heatmap_default_columns_annotations
from .anndata_graphs import ad_add_genes_heatmap_default_rows_annotations
from .anndata_graphs import ad_gene_gene_graph
from .anndata_graphs import ad_genes_heatmap_graph
from .data_sources import fill_block
from .data_sources import fill_boolean
from .data_sources import fill_gene_expression
from .data_sources import fill_global_flow_order
from .data_sources import fill_is_marker
from .data_sources import fill_mean_cells_per_metacell
from .data_sources import fill_mean_total_UMIs_per_cell
from .data_sources import fill_mean_total_UMIs_per_metacell
from .data_sources import fill_n_cells
from .data_sources import fill_n_metacells
from .data_sources import fill_total_UMIs
from .data_sources import fill_type
from .data_sources import fill_umap
from .data_sources import get_is_in_categories_vector
from .data_sources import get_is_in_environment_vector
from .data_sources import get_is_in_neighborhood_vector
from .data_sources import get_skeleton_gene_indices
from .data_sources import get_top_marker_gene_indices
from .heatmap_graphs import add_genes_heatmap_default_columns_annotations
from .heatmap_graphs import add_genes_heatmap_default_rows_annotations
from .heatmap_graphs import genes_heatmap_graph
from .julia_import import jl
from .scatter_graphs import gene_gene_graph

__all__ = ["run_workload"]


def run_workload() -> None:
    """
    Make the calls of a typical session, making each graph from a small fixed ``daf`` chain and ``AnnData``, filling it
    by each data source a notebook offers, and showing it (for internal use).
    """
    with TemporaryDirectory() as directory:
        _run_daf_workload(directory)
        _run_anndata_workload(directory)


def _run_daf_workload(directory: str) -> None:
    # The same calls as ``MetacellsGraphs``' own workload, on the same repositories, opened as a notebook opens them.
    jl.MetacellsGraphs.Workloads.write_workload_daf(directory)
    daf = dp.complete_daf(os.path.join(directory, "top"))

    for axis in ("metacell", "block"):
        points = gene_gene_graph(daf, axis=axis, x_gene="A", y_gene="B")
        fill_type(points.points_colors_vector_fields(), daf, axis=axis)
        fill_n_cells(points.points_sizes_vector_fields(), daf, axis=axis)
        fill_total_UMIs(points.points_entities(), daf, axis=axis)
        fill_mean_total_UMIs_per_cell(points.points_entities(), daf, axis=axis)
        fill_gene_expression(points.points_entities(), daf, axis=axis, gene="C")
        fill_global_flow_order(points.points_entities(), daf, axis=axis)
        _ = points.figure

        for heatmap in (
            genes_heatmap_graph(daf, axis=axis),
            genes_heatmap_graph(
                daf,
                axis=axis,
                genes=get_top_marker_gene_indices(daf, markers_count=3),
                entries=get_is_in_neighborhood_vector(daf, axis=axis, block="B1")
                | get_is_in_categories_vector(daf, axis=axis, property="type", categories=["T2"]),
            ),
        ):
            add_genes_heatmap_default_rows_annotations(heatmap, daf)
            add_genes_heatmap_default_columns_annotations(heatmap, daf, axis=axis)
            fill_is_marker(heatmap.rows_annotations_colors_vector_fields(heatmap.add_rows_annotation()), daf)
            _ = heatmap.figure

    points = gene_gene_graph(daf, axis="metacell", x_gene="A", y_gene="B")
    fill_block(points.points_colors_vector_fields(), daf, axis="metacell")
    fill_boolean(points.points_entities(), daf, axis="metacell", property="is_doublet")
    fill_umap(points.points_entities(), daf, axis="metacell", coordinate="x")
    get_is_in_environment_vector(daf, axis="metacell", block="B1")
    get_skeleton_gene_indices(daf)
    _ = points.figure

    points = gene_gene_graph(daf, axis="block", x_gene="A", y_gene="B")
    fill_n_metacells(points.points_sizes_vector_fields(), daf, axis="block")
    fill_mean_cells_per_metacell(points.points_entities(), daf, axis="block")
    fill_mean_total_UMIs_per_metacell(points.points_entities(), daf, axis="block")
    _ = points.figure


def _run_anndata_workload(directory: str) -> None:
    # A metacells ``AnnData`` in the Python metacells layout, with no marker rank, so the heatmap computes it.
    adata = AnnData(
        X=np.array(
            [[0.1, 0.4, 0.0, 0.2], [0.2, 0.3, 0.1, 0.2], [0.3, 0.2, 0.0, 0.2], [0.4, 0.1, 0.1, 0.2]], dtype="float32"
        ),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(["T1", "T1", "T2", None]),
                "total_umis": np.array([100, 200, 300, 400], dtype="float32"),
                "x": np.array([0, 1, 2, 3], dtype="float32"),
                "y": np.array([3, 2, 1, 0], dtype="float32"),
            },
            index=["M1", "M2", "M3", "M4"],
        ),
        var=pd.DataFrame(
            {
                "lateral_gene": [False, True, False, False],
                "excluded_gene": [False, False, False, True],
                "marker_gene": [True, True, True, False],
            },
            index=["A", "B", "C", "D"],
        ),
    )
    type_colors_csv = os.path.join(directory, "type_colors.csv")
    with open(type_colors_csv, "w", encoding="utf8") as file:
        file.write("type,color\nT1,red\nT2,blue\n")

    points = ad_gene_gene_graph(adata, x_gene="A", y_gene="B")
    ad_fill_type(
        points.points_colors_vector_fields(), adata, type_property="cell_type", type_colors_csv=type_colors_csv
    )
    ad_fill_total_UMIs(points.points_sizes_vector_fields(), adata)
    ad_fill_gene_expression(points.points_entities(), adata, gene="C")
    ad_fill_umap(points.points_entities(), adata, coordinate="x")
    _ = points.figure

    for heatmap in (
        ad_genes_heatmap_graph(adata),
        ad_genes_heatmap_graph(adata, genes=ad_get_top_marker_gene_indices(adata, markers_count=2)),
    ):
        ad_add_genes_heatmap_default_rows_annotations(heatmap, adata)
        ad_add_genes_heatmap_default_columns_annotations(
            heatmap, adata, type_property="cell_type", type_colors_csv=type_colors_csv
        )
        annotation = heatmap.rows_annotations_colors_vector_fields(heatmap.add_rows_annotation())
        ad_fill_is_marker(annotation, adata)
        annotation = heatmap.rows_annotations_colors_vector_fields(heatmap.add_rows_annotation())
        ad_fill_is_excluded(annotation, adata)
        _ = heatmap.figure
