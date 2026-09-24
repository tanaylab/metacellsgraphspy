"""
Data sources for the graphs of a metacells repository. See the Julia
`documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html>`__ for details.

A data source fetches data from a ``Daf`` repository (or a ``pandas`` frame) and writes it into the data source views
of a ``somegraphspy`` graph. The ``get_`` functions fetch and return the data, the ``put_`` functions write data you
already have into the ``sinks``, and the ``fill_`` functions do both. The ``sinks`` are any graph struct (a view, or
one of the data or configuration structs a view holds), or a sequence of them.

Indices are 1-based, as in Julia: the ``entries`` and ``genes`` given by index, the ``indices`` of the low level
queries, and the gene indices returned by the gene selection functions.
"""

# A function per Julia data source, so the module is as long as the Julia one.
# pylint: disable=too-many-lines

from typing import Mapping
from typing import Optional
from typing import Sequence
from typing import Union

import numpy as np
import pandas as pd
from dafpy import DafReader
from somegraphspy import BoolsVector
from somegraphspy import IntegersVector
from somegraphspy import MatrixDataSinks
from somegraphspy import NumbersMatrix
from somegraphspy import NumbersVector
from somegraphspy import Sinks
from somegraphspy import StringsVector
from somegraphspy import VectorDataSinks
from somegraphspy.julia_import import _to_julia

from .julia_import import _from_julia_array
from .julia_import import _given
from .julia_import import _to_julia_array
from .julia_import import _to_julia_frame
from .julia_import import jl

__all__ = [
    "EMPTY_TYPE_COLOR",
    "Entries",
    "GENE_FRACTION_REGULARIZATION_FOR_GRAPHS",
    "MAX_FOLD_FOR_GRAPHS",
    "Via",
    "fill_axes_matrix_data",
    "fill_axes_names_data",
    "fill_axis_names_data",
    "fill_axis_vector_data",
    "fill_block",
    "fill_boolean_annotation",
    "fill_column_boolean_annotation",
    "fill_column_names_data",
    "fill_column_vector_data",
    "fill_gene_correlation",
    "fill_gene_correlation_change",
    "fill_gene_expression",
    "fill_genes_expression_matrix",
    "fill_genes_fold_matrix",
    "fill_global_flow_order",
    "fill_mean_cells_per_metacell",
    "fill_mean_total_UMIs_per_cell",
    "fill_mean_total_UMIs_per_metacell",
    "fill_module_regulators_hovers",
    "fill_n_cells",
    "fill_n_metacells",
    "fill_total_UMIs",
    "fill_type",
    "fill_umap",
    "get_axes_matrix",
    "get_axis_entries_vector",
    "get_axis_vector",
    "get_block_vector",
    "get_boolean_annotation_vector",
    "get_column_vector",
    "get_gene_correlation_change_vector",
    "get_gene_correlation_vector",
    "get_gene_expression_vector",
    "get_genes_expression_matrix",
    "get_genes_fold_matrix",
    "get_global_flow_order_vector",
    "get_matrix_query",
    "get_mean_cells_per_metacell_vector",
    "get_mean_total_UMIs_per_cell_vector",
    "get_mean_total_UMIs_per_metacell_vector",
    "get_module_regulators_hovers",
    "get_n_cells_vector",
    "get_n_metacells_vector",
    "get_skeleton_gene_indices",
    "get_top_marker_gene_indices",
    "get_total_UMIs_vector",
    "get_type_colors",
    "get_type_vector",
    "get_umap_vector",
    "get_vector_query",
    "put_boolean_annotation_configuration",
    "put_count_configuration",
    "put_gene_correlation_change_configuration",
    "put_genes_expression_configuration",
    "put_genes_fold_configuration",
    "put_matrix_data",
    "put_matrix_names_data",
    "put_type_configuration",
    "put_umap_configuration",
    "put_umap_data",
    "put_vector_data",
    "put_vector_mask_data",
    "put_vector_names_data",
]

#: Which entries of an axis a data source uses: their names, or their (1-based) indices.
Entries = Union[StringsVector, IntegersVector]

#: A chain of axes to walk before asking for a property (e.g., ``["metacell"]`` to ask for the block of the metacell of
#: each cell).
Via = Sequence[str]

#: The regularization added to a gene fraction before taking its log. See the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.GENE_FRACTION_REGULARIZATION_FOR_GRAPHS>`__
#: for details.
GENE_FRACTION_REGULARIZATION_FOR_GRAPHS: float = float(jl.MetacellsGraphs.GENE_FRACTION_REGULARIZATION_FOR_GRAPHS)

#: The extreme of the color scale of the log-base-2 fold factors between gene expression and median gene expression.
#: See the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.MAX_FOLD_FOR_GRAPHS>`__
#: for details.
MAX_FOLD_FOR_GRAPHS: float = float(jl.MetacellsGraphs.MAX_FOLD_FOR_GRAPHS)

#: The color to give to entities without any type annotation. See the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.EMPTY_TYPE_COLOR>`__
#: for details.
EMPTY_TYPE_COLOR: str = str(jl.MetacellsGraphs.EMPTY_TYPE_COLOR)


# Generic axis names.


def fill_axis_names_data(
    sinks: VectorDataSinks, daf: DafReader, *, axis: str, entries: Optional[Entries] = None
) -> None:
    """
    Name the entities of the ``sinks`` after each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_axis_names_data!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_axis_names_data_b(
        _to_julia(sinks), daf, axis=axis, **_given(entries=_to_julia_array(entries))
    )


def put_vector_names_data(sinks: VectorDataSinks, name_per_entry: StringsVector) -> None:
    """
    Name the entities of the ``sinks`` after the ``name_per_entry``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_vector_names_data!>`__
    for details.
    """
    jl.MetacellsGraphs.put_vector_names_data_b(_to_julia(sinks), _to_julia_array(name_per_entry))


def get_axis_entries_vector(daf: DafReader, *, axis: str, entries: Optional[Entries] = None) -> np.ndarray:
    """
    Get the name for each of the ``entries`` of the ``daf`` ``axis`` (by default, all of them). See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_axis_entries_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_axis_entries_vector(daf, axis=axis, **_given(entries=_to_julia_array(entries)))
    )


# Generic axis vectors.


def fill_axis_vector_data(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: str,
    query_suffix: str,
    via: Optional[Via] = None,
    empty_value: Optional[Union[str, float, bool]] = None,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
) -> None:
    """
    Fill the ``sinks`` with the result of the query ``@ axis query_suffix`` for each of the ``entries`` of the ``daf``
    ``axis``, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_axis_vector_data!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_axis_vector_data_b(
        _to_julia(sinks),
        daf,
        axis=axis,
        query_suffix=query_suffix,
        **_given(via=_to_julia_array(via), empty_value=empty_value, entries=_to_julia_array(entries), title=title),
    )


def put_vector_data(
    sinks: VectorDataSinks,
    value_per_entry: Union[NumbersVector, StringsVector, BoolsVector],
    *,
    title: Optional[str] = None,
) -> None:
    """
    Put a ``value_per_entry`` of an axis into the ``sinks``: as the values of a role, and as a hover line on the
    entities. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_vector_data!>`__
    for details.
    """
    jl.MetacellsGraphs.put_vector_data_b(_to_julia(sinks), _to_julia_array(value_per_entry), **_given(title=title))


def get_axis_vector(
    daf: DafReader,
    *,
    axis: str,
    query_suffix: str,
    via: Optional[Via] = None,
    empty_value: Optional[Union[str, float, bool]] = None,
    entries: Optional[Entries] = None,
) -> np.ndarray:
    """
    Get the result of the query ``@ axis via query_suffix`` for each of the ``entries`` of the ``daf`` ``axis``. See the
    Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_axis_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_axis_vector(
            daf,
            axis=axis,
            query_suffix=query_suffix,
            **_given(via=_to_julia_array(via), empty_value=empty_value, entries=_to_julia_array(entries)),
        )
    )


# Masks.


def put_vector_mask_data(sinks: VectorDataSinks, is_shown_per_entry: BoolsVector) -> None:
    """
    Hide the entities of the ``sinks`` which are not shown by the ``is_shown_per_entry`` mask. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_vector_mask_data!>`__
    for details.
    """
    jl.MetacellsGraphs.put_vector_mask_data_b(_to_julia(sinks), _to_julia_array(is_shown_per_entry))


# Frame columns.


def fill_column_vector_data(
    sinks: VectorDataSinks, frame: pd.DataFrame, *, column: str, title: Optional[str] = None
) -> None:
    """
    Fill the ``sinks`` with the ``column`` of a ``frame``, one value per row. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_column_vector_data!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_column_vector_data_b(
        _to_julia(sinks), _to_julia_frame(frame), column=column, **_given(title=title)
    )


def fill_column_names_data(sinks: VectorDataSinks, frame: pd.DataFrame, *, column: str) -> None:
    """
    Name the entities of the ``sinks`` after the ``column`` of a ``frame``, one name per row. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_column_names_data!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_column_names_data_b(_to_julia(sinks), _to_julia_frame(frame), column=column)


def fill_column_boolean_annotation(
    sinks: VectorDataSinks,
    frame: pd.DataFrame,
    *,
    column: str,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the Boolean ``column`` of a ``frame``, shown as ``true`` in black and ``false`` in light
    grey. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_column_boolean_annotation!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_column_boolean_annotation_b(
        _to_julia(sinks), _to_julia_frame(frame), column=column, **_given(title=title, show_legend=show_legend)
    )


def get_column_vector(frame: pd.DataFrame, *, column: str) -> np.ndarray:
    """
    Get the ``column`` of a ``frame``, one value per row. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_column_vector>`__
    for details.
    """
    return _from_julia_array(jl.MetacellsGraphs.get_column_vector(_to_julia_frame(frame), column=column))


# Module regulators.


def fill_module_regulators_hovers(
    sinks: VectorDataSinks, frame: pd.DataFrame, *, prefix: str, side_name: str, regulators_count: int
) -> None:
    """
    Add a hover per row of a gene report ``frame`` to the ``sinks``, saying which regulators the gene is in a module
    with in the base blocks of one side. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_module_regulators_hovers!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_module_regulators_hovers_b(
        _to_julia(sinks), _to_julia_frame(frame), prefix=prefix, side_name=side_name, regulators_count=regulators_count
    )


def get_module_regulators_hovers(
    frame: pd.DataFrame, *, prefix: str, side_name: str, regulators_count: int
) -> np.ndarray:
    """
    Get a hover per row of a gene report ``frame``, saying how often the gene is in no module at all in the base blocks
    of one side, and which regulators it is most often in a module with there. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_module_regulators_hovers>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_module_regulators_hovers(
            _to_julia_frame(frame), prefix=prefix, side_name=side_name, regulators_count=regulators_count
        )
    )


# Generic axes names and matrices.


def fill_axes_names_data(
    sinks: MatrixDataSinks,
    daf: DafReader,
    *,
    rows_axis: str,
    columns_axis: str,
    row_entries: Optional[Entries] = None,
    column_entries: Optional[Entries] = None,
) -> None:
    """
    Name the rows of the ``sinks`` after the ``row_entries`` of the ``daf`` ``rows_axis``, and their columns after the
    ``column_entries`` of its ``columns_axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_axes_names_data!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_axes_names_data_b(
        _to_julia(sinks),
        daf,
        rows_axis=rows_axis,
        columns_axis=columns_axis,
        **_given(row_entries=_to_julia_array(row_entries), column_entries=_to_julia_array(column_entries)),
    )


def put_matrix_names_data(sinks: MatrixDataSinks, name_per_row: StringsVector, name_per_column: StringsVector) -> None:
    """
    Name the rows and the columns of the ``sinks`` after the ``name_per_row`` and the ``name_per_column``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_matrix_names_data!>`__
    for details.
    """
    jl.MetacellsGraphs.put_matrix_names_data_b(
        _to_julia(sinks), _to_julia_array(name_per_row), _to_julia_array(name_per_column)
    )


def fill_axes_matrix_data(
    sinks: MatrixDataSinks,
    daf: DafReader,
    *,
    rows_axis: str,
    columns_axis: str,
    query_suffix: str,
    row_entries: Optional[Entries] = None,
    column_entries: Optional[Entries] = None,
    title: Optional[str] = None,
) -> None:
    """
    Fill the ``sinks`` with the result of the query ``@ rows_axis @ columns_axis query_suffix``, and name the rows and
    the columns. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_axes_matrix_data!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_axes_matrix_data_b(
        _to_julia(sinks),
        daf,
        rows_axis=rows_axis,
        columns_axis=columns_axis,
        query_suffix=query_suffix,
        **_given(row_entries=_to_julia_array(row_entries), column_entries=_to_julia_array(column_entries), title=title),
    )


def put_matrix_data(
    sinks: MatrixDataSinks, value_per_row_per_column: NumbersMatrix, *, title: Optional[str] = None
) -> None:
    """
    Put a ``value_per_row_per_column`` of two axes into the ``sinks``: as the values of the entries, and as a hover line
    on each entry. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_matrix_data!>`__
    for details.
    """
    jl.MetacellsGraphs.put_matrix_data_b(
        _to_julia(sinks), _to_julia_array(value_per_row_per_column), **_given(title=title)
    )


def get_axes_matrix(
    daf: DafReader,
    *,
    rows_axis: str,
    columns_axis: str,
    query_suffix: str,
    row_entries: Optional[Entries] = None,
    column_entries: Optional[Entries] = None,
) -> np.ndarray:
    """
    Get the result of the query ``@ rows_axis @ columns_axis query_suffix`` for each of the ``row_entries`` and
    ``column_entries`` of the two ``daf`` axes. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_axes_matrix>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_axes_matrix(
            daf,
            rows_axis=rows_axis,
            columns_axis=columns_axis,
            query_suffix=query_suffix,
            **_given(row_entries=_to_julia_array(row_entries), column_entries=_to_julia_array(column_entries)),
        )
    )


# Gene expression.


def put_genes_expression_configuration(
    sinks: Sinks,
    *,
    gene_fraction_regularization: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Show gene expression (linear fraction) on a log-base-2 scale, regularized by ``gene_fraction_regularization``,
    named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_genes_expression_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_genes_expression_configuration_b(
        _to_julia(sinks),
        **_given(gene_fraction_regularization=gene_fraction_regularization, title=title, show_legend=show_legend),
    )


def fill_gene_expression(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    gene: str,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    gene_fraction_regularization: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the ``gene`` expression level (linear fraction) per each of the ``entries`` of some ``daf``
    ``axis``, shown in log base 2, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_gene_expression!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_gene_expression_b(
        _to_julia(sinks),
        daf,
        gene=gene,
        **_given(
            axis=axis,
            entries=_to_julia_array(entries),
            gene_fraction_regularization=gene_fraction_regularization,
            title=title,
            show_legend=show_legend,
        ),
    )


def get_gene_expression_vector(
    daf: DafReader, *, gene: str, axis: Optional[str] = None, entries: Optional[Entries] = None
) -> np.ndarray:
    """
    Get the vector of the ``gene`` expression level (linear fraction) per each of the ``entries`` of some ``daf``
    ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_gene_expression_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_gene_expression_vector(
            daf, gene=gene, **_given(axis=axis, entries=_to_julia_array(entries))
        )
    )


def fill_genes_expression_matrix(
    sinks: MatrixDataSinks,
    daf: DafReader,
    *,
    genes: Optional[Entries] = None,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    gene_fraction_regularization: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the ``genes`` expression level (linear fraction) per each of the ``entries`` of some ``daf``
    ``axis``, shown in log base 2, and name the rows and the columns. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_genes_expression_matrix!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_genes_expression_matrix_b(
        _to_julia(sinks),
        daf,
        **_given(
            genes=_to_julia_array(genes),
            axis=axis,
            entries=_to_julia_array(entries),
            gene_fraction_regularization=gene_fraction_regularization,
            title=title,
            show_legend=show_legend,
        ),
    )


def get_genes_expression_matrix(
    daf: DafReader, *, genes: Optional[Entries] = None, axis: Optional[str] = None, entries: Optional[Entries] = None
) -> np.ndarray:
    """
    Get the matrix of the ``genes`` expression level (linear fraction) per each of the ``entries`` of some ``daf``
    ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_genes_expression_matrix>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_genes_expression_matrix(
            daf, **_given(genes=_to_julia_array(genes), axis=axis, entries=_to_julia_array(entries))
        )
    )


# Gene fold matrices.


def fill_genes_fold_matrix(
    sinks: MatrixDataSinks,
    daf: DafReader,
    *,
    genes: Optional[Entries] = None,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    max_fold: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the ``genes`` fold (log base-2 minus the median) per each of the ``entries`` of some ``daf``
    ``axis``, shown in the range -``max_fold`` (blue) to 0 (white) to +``max_fold`` (red), and name the rows and the
    columns. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_genes_fold_matrix!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_genes_fold_matrix_b(
        _to_julia(sinks),
        daf,
        **_given(
            genes=_to_julia_array(genes),
            axis=axis,
            entries=_to_julia_array(entries),
            max_fold=max_fold,
            title=title,
            show_legend=show_legend,
        ),
    )


def put_genes_fold_configuration(
    sinks: Sinks, *, max_fold: Optional[float] = None, title: Optional[str] = None, show_legend: Optional[bool] = None
) -> None:
    """
    Show the genes fold (log base-2 minus the median) in the range -``max_fold`` (blue) to 0 (white) to +``max_fold``
    (red), named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_genes_fold_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_genes_fold_configuration_b(
        _to_julia(sinks), **_given(max_fold=max_fold, title=title, show_legend=show_legend)
    )


def get_genes_fold_matrix(
    daf: DafReader, *, genes: Optional[Entries] = None, axis: Optional[str] = None, entries: Optional[Entries] = None
) -> np.ndarray:
    """
    Get the matrix of the ``genes`` fold (log base-2 minus the median) per each of the ``entries`` of some ``daf``
    ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_genes_fold_matrix>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_genes_fold_matrix(
            daf, **_given(genes=_to_julia_array(genes), axis=axis, entries=_to_julia_array(entries))
        )
    )


# Gene correlations.


def fill_gene_correlation(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    gene: str,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
) -> None:
    """
    Fill the ``sinks`` with the correlation of the ``gene`` per each of the ``entries`` of some ``daf`` ``axis``, and
    name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_gene_correlation!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_gene_correlation_b(
        _to_julia(sinks), daf, gene=gene, **_given(axis=axis, entries=_to_julia_array(entries), title=title)
    )


def get_gene_correlation_vector(
    daf: DafReader, *, gene: str, axis: Optional[str] = None, entries: Optional[Entries] = None
) -> np.ndarray:
    """
    Get the correlation of the ``gene`` between the cells of the neighborhood of each of the ``entries`` of some
    ``daf`` ``axis`` and their punctuated metacells. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_gene_correlation_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_gene_correlation_vector(
            daf, gene=gene, **_given(axis=axis, entries=_to_julia_array(entries))
        )
    )


def fill_gene_correlation_change(
    sinks: VectorDataSinks,
    daf: DafReader,
    base_daf: DafReader,
    *,
    gene: str,
    axis: Optional[str] = None,
    base_axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with how much the correlation of the ``gene`` changed between the ``base_daf`` and the ``daf``,
    per each of the ``entries``, and name the entities after them. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_gene_correlation_change!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_gene_correlation_change_b(
        _to_julia(sinks),
        daf,
        base_daf,
        gene=gene,
        **_given(
            axis=axis, base_axis=base_axis, entries=_to_julia_array(entries), title=title, show_legend=show_legend
        ),
    )


def put_gene_correlation_change_configuration(
    sinks: Sinks, *, title: Optional[str] = None, show_legend: Optional[bool] = None
) -> None:
    """
    Name the correlation change in the configuration, and include it in the legend if ``show_legend`` where it is the
    colors. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_gene_correlation_change_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_gene_correlation_change_configuration_b(
        _to_julia(sinks), **_given(title=title, show_legend=show_legend)
    )


def get_gene_correlation_change_vector(
    daf: DafReader,
    base_daf: DafReader,
    *,
    gene: str,
    axis: Optional[str] = None,
    base_axis: Optional[str] = None,
    entries: Optional[Entries] = None,
) -> np.ndarray:
    """
    Get how much the correlation of the ``gene`` changed between the ``base_daf`` and the ``daf``, per each of the
    ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_gene_correlation_change_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_gene_correlation_change_vector(
            daf, base_daf, gene=gene, **_given(axis=axis, base_axis=base_axis, entries=_to_julia_array(entries))
        )
    )


# UMAP coordinates.


def fill_umap(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    coordinate: str,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
) -> None:
    """
    Fill the ``sinks`` with the UMAP ``coordinate`` per each of the ``entries`` of some ``daf`` ``axis``, shown without
    ticks or a grid, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_umap!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_umap_b(
        _to_julia(sinks),
        daf,
        coordinate=coordinate,
        **_given(axis=axis, entries=_to_julia_array(entries), title=title),
    )


def put_umap_data(sinks: VectorDataSinks, coordinate_per_entry: NumbersVector, *, title: Optional[str] = None) -> None:
    """
    Put a UMAP ``coordinate_per_entry`` into the ``sinks``, as the values of a role. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_umap_data!>`__
    for details.
    """
    jl.MetacellsGraphs.put_umap_data_b(_to_julia(sinks), _to_julia_array(coordinate_per_entry), **_given(title=title))


def put_umap_configuration(sinks: Sinks, *, title: Optional[str] = None) -> None:
    """
    Show UMAP coordinates without ticks or a grid, named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_umap_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_umap_configuration_b(_to_julia(sinks), **_given(title=title))


def get_umap_vector(
    daf: DafReader, *, coordinate: str, axis: Optional[str] = None, entries: Optional[Entries] = None
) -> np.ndarray:
    """
    Get the vector of the UMAP ``coordinate`` per each of the ``entries`` of some ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_umap_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_umap_vector(
            daf, coordinate=coordinate, **_given(axis=axis, entries=_to_julia_array(entries))
        )
    )


# Gene selection.


def get_top_marker_gene_indices(daf: DafReader, *, markers_count: int) -> np.ndarray:
    """
    Get the (1-based) indices of the ``markers_count`` best marker genes of a ``daf`` repository, to pass as the
    ``genes`` of another data source. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_top_marker_gene_indices>`__
    for details.
    """
    return _from_julia_array(jl.MetacellsGraphs.get_top_marker_gene_indices(daf, markers_count=markers_count))


def get_skeleton_gene_indices(daf: DafReader) -> np.ndarray:
    """
    Get the (1-based) indices of the skeleton genes of a ``daf`` repository, to pass as the ``genes`` of another data
    source. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_skeleton_gene_indices>`__
    for details.
    """
    return _from_julia_array(jl.MetacellsGraphs.get_skeleton_gene_indices(daf))


# Counts.


def put_count_configuration(sinks: Sinks, *, title: Optional[str] = None, show_legend: Optional[bool] = None) -> None:
    """
    Show a count on a log-base-2 scale, in the ``YlOrRd`` color scale, named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_count_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_count_configuration_b(_to_julia(sinks), **_given(title=title, show_legend=show_legend))


def fill_total_UMIs(  # pylint: disable=invalid-name
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the total UMIs of each of the ``entries`` of the ``daf`` ``axis``, shown as a count, and
    name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_total_UMIs!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_total_UMIs_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis,
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            empty_value=empty_value,
            title=title,
            show_legend=show_legend,
        ),
    )


def get_total_UMIs_vector(  # pylint: disable=invalid-name
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
) -> np.ndarray:
    """
    Get the total UMIs of each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_total_UMIs_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_total_UMIs_vector(
            daf,
            **_given(axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value),
        )
    )


def fill_n_cells(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the number of cells of each of the ``entries`` of the ``daf`` ``axis``, shown as a count,
    and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_n_cells!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_n_cells_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis,
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            empty_value=empty_value,
            title=title,
            show_legend=show_legend,
        ),
    )


def get_n_cells_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
) -> np.ndarray:
    """
    Get the number of cells of each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_n_cells_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_n_cells_vector(
            daf,
            **_given(axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value),
        )
    )


def fill_n_metacells(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the number of metacells of each of the ``entries`` of the ``daf`` ``axis`` (by default, the
    blocks), shown as a count, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_n_metacells!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_n_metacells_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis,
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            empty_value=empty_value,
            title=title,
            show_legend=show_legend,
        ),
    )


def get_n_metacells_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
) -> np.ndarray:
    """
    Get the number of metacells of each of the ``entries`` of the ``daf`` ``axis`` (by default, the blocks). See the
    Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_n_metacells_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_n_metacells_vector(
            daf,
            **_given(axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value),
        )
    )


def fill_mean_cells_per_metacell(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the mean number of cells per metacell of each of the ``entries`` of the ``daf`` ``axis``
    (by default, the blocks), shown as a count, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_mean_cells_per_metacell!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_mean_cells_per_metacell_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), title=title, show_legend=show_legend
        ),
    )


def get_mean_cells_per_metacell_vector(
    daf: DafReader, *, axis: Optional[str] = None, entries: Optional[Entries] = None, via: Optional[Via] = None
) -> np.ndarray:
    """
    Get the mean number of cells per metacell of each of the ``entries`` of the ``daf`` ``axis`` (by default, the
    blocks). See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_mean_cells_per_metacell_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_mean_cells_per_metacell_vector(
            daf, **_given(axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via))
        )
    )


def fill_mean_total_UMIs_per_metacell(  # pylint: disable=invalid-name
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the mean total UMIs per metacell of each of the ``entries`` of the ``daf`` ``axis`` (by
    default, the blocks), shown as a count, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_mean_total_UMIs_per_metacell!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_mean_total_UMIs_per_metacell_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), title=title, show_legend=show_legend
        ),
    )


def get_mean_total_UMIs_per_metacell_vector(  # pylint: disable=invalid-name
    daf: DafReader, *, axis: Optional[str] = None, entries: Optional[Entries] = None, via: Optional[Via] = None
) -> np.ndarray:
    """
    Get the mean total UMIs per metacell of each of the ``entries`` of the ``daf`` ``axis`` (by default, the blocks).
    See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_mean_total_UMIs_per_metacell_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_mean_total_UMIs_per_metacell_vector(
            daf, **_given(axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via))
        )
    )


def fill_mean_total_UMIs_per_cell(  # pylint: disable=invalid-name
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the mean total UMIs per cell of each of the ``entries`` of the ``daf`` ``axis``, shown as a
    count, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_mean_total_UMIs_per_cell!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_mean_total_UMIs_per_cell_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), title=title, show_legend=show_legend
        ),
    )


def get_mean_total_UMIs_per_cell_vector(  # pylint: disable=invalid-name
    daf: DafReader, *, axis: Optional[str] = None, entries: Optional[Entries] = None, via: Optional[Via] = None
) -> np.ndarray:
    """
    Get the mean total UMIs per cell of each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_mean_total_UMIs_per_cell_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_mean_total_UMIs_per_cell_vector(
            daf, **_given(axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via))
        )
    )


# Blocks.


def fill_block(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[str] = None,
    title: Optional[str] = None,
) -> None:
    """
    Fill the ``sinks`` with the block of each of the ``entries`` of the ``daf`` ``axis``, and name the entities after
    the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_block!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_block_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value, title=title
        ),
    )


def get_block_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[str] = None,
) -> np.ndarray:
    """
    Get the block of each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_block_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_block_vector(
            daf,
            **_given(axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value),
        )
    )


# Global flow order.


def fill_global_flow_order(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    type_property: Optional[str] = None,
    empty_value: Optional[float] = None,
    title: Optional[str] = None,
) -> None:
    """
    Fill the ``sinks`` with the global flow order of the type of each of the ``entries`` of the ``daf`` ``axis``, and
    name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_global_flow_order!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_global_flow_order_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis,
            entries=_to_julia_array(entries),
            type_property=type_property,
            empty_value=empty_value,
            title=title,
        ),
    )


def get_global_flow_order_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    type_property: Optional[str] = None,
    empty_value: Optional[float] = None,
) -> np.ndarray:
    """
    Get the global flow order of the type of each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_global_flow_order_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_global_flow_order_vector(
            daf,
            **_given(axis=axis, entries=_to_julia_array(entries), type_property=type_property, empty_value=empty_value),
        )
    )


# Types.


def fill_type(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    type_property: Optional[str] = None,
    type_axis: Optional[str] = None,
    via: Optional[Via] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
    empty_type_color: Optional[str] = None,
) -> None:
    """
    Fill the ``sinks`` with the ``type_property`` of each of the ``entries`` of the ``daf`` ``axis``, colored by the
    ``color`` of the ``type_axis``, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_type!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_type_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis,
            entries=_to_julia_array(entries),
            type_property=type_property,
            type_axis=type_axis,
            via=_to_julia_array(via),
            title=title,
            show_legend=show_legend,
            empty_type_color=empty_type_color,
        ),
    )


def put_type_configuration(
    sinks: Sinks, color_per_type: Mapping[str, str], *, title: Optional[str] = None, show_legend: Optional[bool] = None
) -> None:
    """
    Show the types of the ``sinks`` using the ``color_per_type`` palette, named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_type_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_type_configuration_b(
        _to_julia(sinks), _to_julia(dict(color_per_type)), **_given(title=title, show_legend=show_legend)
    )


def get_type_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    type_property: Optional[str] = None,
    via: Optional[Via] = None,
) -> np.ndarray:
    """
    Get the ``type_property`` for each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_type_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_type_vector(
            daf,
            **_given(
                axis=axis, entries=_to_julia_array(entries), type_property=type_property, via=_to_julia_array(via)
            ),
        )
    )


def get_type_colors(
    daf: DafReader, *, type_axis: Optional[str] = None, empty_type_color: Optional[str] = None
) -> Mapping[str, str]:
    """
    Get the palette mapping each entry of the ``daf`` ``type_axis`` to its ``color``, with an additional
    ``empty_type_color`` for entries with no type (empty string). See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_type_colors>`__
    for details.
    """
    colors = jl.MetacellsGraphs.get_type_colors(daf, **_given(type_axis=type_axis, empty_type_color=empty_type_color))
    return {str(key): str(jl.getindex(colors, key)) for key in jl.keys(colors)}


# Boolean annotations.


def fill_boolean_annotation(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the Boolean mask ``property`` of each of the ``entries`` of the ``daf`` ``axis``, shown as
    ``true`` in black and ``false`` in light grey, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_boolean_annotation!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_boolean_annotation_b(
        _to_julia(sinks),
        daf,
        axis=axis,
        property=property,
        **_given(entries=_to_julia_array(entries), title=title, show_legend=show_legend),
    )


def put_boolean_annotation_configuration(
    sinks: Sinks, *, title: Optional[str] = None, show_legend: Optional[bool] = None
) -> None:
    """
    Show Boolean *string* annotations: ``true`` in black and ``false`` in light grey, named by the ``title``. See the
    Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_boolean_annotation_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_boolean_annotation_configuration_b(
        _to_julia(sinks), **_given(title=title, show_legend=show_legend)
    )


def get_boolean_annotation_vector(
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
) -> np.ndarray:
    """
    Get the Boolean mask ``property`` of each of the ``entries`` of the ``daf`` ``axis``, as the *strings* ``"true"``
    and ``"false"``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_boolean_annotation_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_boolean_annotation_vector(
            daf, axis=axis, property=property, **_given(entries=_to_julia_array(entries))
        )
    )


# Low level queries.


def get_vector_query(daf: DafReader, *, query: str, indices: Optional[Sequence[int]] = None) -> np.ndarray:
    """
    Execute a vector ``query`` on a ``daf`` repository and optionally access specific (1-based) ``indices`` in it. See
    the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_vector_query>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_vector_query(daf, query=query, **_given(indices=_to_julia_array(indices)))
    )


def get_matrix_query(
    daf: DafReader,
    *,
    query: str,
    row_indices: Optional[Sequence[int]] = None,
    column_indices: Optional[Sequence[int]] = None,
) -> np.ndarray:
    """
    Execute a matrix ``query`` on a ``daf`` repository and optionally access specific (1-based) ``row_indices`` and/or
    ``column_indices`` in it. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_matrix_query>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_matrix_query(
            daf,
            query=query,
            **_given(row_indices=_to_julia_array(row_indices), column_indices=_to_julia_array(column_indices)),
        )
    )
