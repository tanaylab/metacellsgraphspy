"""
Data sources for the graphs of a metacells repository. See the Julia
`documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html>`__ for details.

A data source fetches data from a ``Daf`` repository (or a ``pandas`` frame) and writes it into the data source views
of a ``somegraphspy`` graph. The ``get_`` functions fetch and return the data, the ``put_`` functions write data you
already have into the ``sinks``, and the ``fill_`` functions do both. The ``sinks`` are any graph struct (a view, or
one of the data or configuration structs a view holds), or a sequence of them.

A ``fill_`` function which isn't given its ``entries`` (or ``genes``) fills the entries its ``sinks`` are already named
after, or else the whole axis. So when the first data source of a graph picks some entries, the rest follow it without
being told which they are.

Indices are 1-based, as in Julia: the ``entries`` and ``genes`` given by index, the ``indices`` of the low level
queries, and the gene indices returned by the gene selection functions.

The ``entries`` and ``genes`` can also be given by a Boolean mask over the axis. In particular, any Boolean get can serve
as the mask (e.g., ``entries=get_is_marker_vector(daf)``).
"""

# A function per Julia data source, so the module is as long as the Julia one, and its wrappers look alike.
# pylint: disable=too-many-lines,duplicate-code

from typing import Any
from typing import Callable
from typing import Mapping
from typing import Optional
from typing import Sequence
from typing import Union

import numpy as np
import pandas as pd
import somegraphspy
from dafpy import DafReader
from somegraphspy import BoolsVector
from somegraphspy import HoverFormat
from somegraphspy import IntegersVector
from somegraphspy import LogBase
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
    "CORRELATION_HOVER_FORMAT_FOR_GRAPHS",
    "EMPTY_CATEGORY_COLOR",
    "EXPRESSION_LEVEL_HOVER_FORMAT_FOR_GRAPHS",
    "Entries",
    "FOLD_HOVER_FORMAT_FOR_GRAPHS",
    "GENE_FRACTION_REGULARIZATION_FOR_GRAPHS",
    "MAX_FOLD_FOR_GRAPHS",
    "MEAN_SIZE_HOVER_FORMAT_FOR_GRAPHS",
    "SIZE_HOVER_FORMAT_FOR_GRAPHS",
    "TOP_MARKERS_COUNT_FOR_GRAPHS",
    "Via",
    "fill_axes_matrix_data",
    "fill_axes_names_data",
    "fill_axis_names_data",
    "fill_axis_vector_data",
    "fill_block",
    "fill_boolean",
    "fill_category",
    "fill_column_boolean_annotation",
    "fill_column_names_data",
    "fill_column_vector_data",
    "fill_correlation",
    "fill_correlation_between_cells_and_projected_metacells",
    "fill_correlation_between_cells_and_projected_punctuated_metacells",
    "fill_correlation_between_cells_and_punctuated_metacells",
    "fill_gene_correlation",
    "fill_gene_correlation_change",
    "fill_gene_expression",
    "fill_genes_expression_matrix",
    "fill_genes_fold_matrix",
    "fill_global_flow_order",
    "fill_is_correlated_with_skeleton",
    "fill_is_excluded",
    "fill_is_forbidden",
    "fill_is_lateral",
    "fill_is_marker",
    "fill_is_mitochondrial",
    "fill_is_regulator",
    "fill_is_ribosomal",
    "fill_is_skeleton",
    "fill_is_transcription_factor",
    "fill_mean_cells_per_metacell",
    "fill_mean_total_UMIs_per_cell",
    "fill_mean_total_UMIs_per_metacell",
    "fill_module_regulators_hovers",
    "fill_n_cells",
    "fill_n_metacells",
    "fill_size",
    "fill_total_UMIs",
    "fill_type",
    "fill_umap",
    "function_default",
    "get_axes_matrix",
    "get_axis_entries_vector",
    "get_axis_vector",
    "get_block_colors",
    "get_block_vector",
    "get_boolean_vector",
    "get_category_colors",
    "get_category_vector",
    "get_column_vector",
    "get_correlation_between_cells_and_projected_metacells_vector",
    "get_correlation_between_cells_and_projected_punctuated_metacells_vector",
    "get_correlation_between_cells_and_punctuated_metacells_vector",
    "get_correlation_vector",
    "get_gene_correlation_change_vector",
    "get_gene_correlation_vector",
    "get_gene_expression_vector",
    "get_gene_is_correlated_vector",
    "get_genes_expression_matrix",
    "get_genes_fold_matrix",
    "get_global_flow_order_vector",
    "get_is_correlated_with_skeleton_vector",
    "get_is_excluded_vector",
    "get_is_forbidden_vector",
    "get_is_in_categories_vector",
    "get_is_in_environment_vector",
    "get_is_in_neighborhood_vector",
    "get_is_lateral_vector",
    "get_is_marker_vector",
    "get_is_mitochondrial_vector",
    "get_is_regulator_vector",
    "get_is_ribosomal_vector",
    "get_is_skeleton_vector",
    "get_is_transcription_factor_vector",
    "get_matrix_query",
    "get_mean_cells_per_metacell_vector",
    "get_mean_total_UMIs_per_cell_vector",
    "get_mean_total_UMIs_per_metacell_vector",
    "get_module_regulators_hovers",
    "get_n_cells_vector",
    "get_n_metacells_vector",
    "get_size_vector",
    "get_skeleton_gene_indices",
    "get_top_marker_gene_indices",
    "get_total_UMIs_vector",
    "get_type_colors",
    "get_type_vector",
    "get_umap_vector",
    "get_vector_query",
    "put_boolean_configuration",
    "put_boolean_values",
    "put_category_configuration",
    "put_category_values",
    "put_coordinates_configuration",
    "put_correlation_configuration",
    "put_correlation_values",
    "put_expression_level_configuration",
    "put_fold_configuration",
    "put_size_configuration",
    "put_size_values",
    "put_umap_data",
]

#: Which entries of an axis a data source uses: their names, their (1-based) indices, or a Boolean mask over the axis.
Entries = Union[StringsVector, IntegersVector, BoolsVector]

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

#: The number of the strongest marker genes a graph shows when it isn't told which genes to show. See the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.TOP_MARKERS_COUNT_FOR_GRAPHS>`__
#: for details.
TOP_MARKERS_COUNT_FOR_GRAPHS: int = int(jl.MetacellsGraphs.TOP_MARKERS_COUNT_FOR_GRAPHS)

#: The ``hover_format`` of a size (a count): a whole number. See the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.SIZE_HOVER_FORMAT_FOR_GRAPHS>`__
#: for details.
SIZE_HOVER_FORMAT_FOR_GRAPHS: HoverFormat = HoverFormat(jl.MetacellsGraphs.SIZE_HOVER_FORMAT_FOR_GRAPHS)

#: The ``hover_format`` of a mean of sizes (e.g., the mean cells per metacell): one digit after the decimal point. See
#: the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.MEAN_SIZE_HOVER_FORMAT_FOR_GRAPHS>`__
#: for details.
MEAN_SIZE_HOVER_FORMAT_FOR_GRAPHS: HoverFormat = HoverFormat(jl.MetacellsGraphs.MEAN_SIZE_HOVER_FORMAT_FOR_GRAPHS)

#: The ``hover_format`` of an expression level: 3 significant digits. See the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.EXPRESSION_LEVEL_HOVER_FORMAT_FOR_GRAPHS>`__
#: for details.
EXPRESSION_LEVEL_HOVER_FORMAT_FOR_GRAPHS: HoverFormat = HoverFormat(
    jl.MetacellsGraphs.EXPRESSION_LEVEL_HOVER_FORMAT_FOR_GRAPHS
)

#: The ``hover_format`` of a fold factor: two digits after the decimal point. See the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.FOLD_HOVER_FORMAT_FOR_GRAPHS>`__
#: for details.
FOLD_HOVER_FORMAT_FOR_GRAPHS: HoverFormat = HoverFormat(jl.MetacellsGraphs.FOLD_HOVER_FORMAT_FOR_GRAPHS)

#: The ``hover_format`` of a correlation: three digits after the decimal point, since small correlations are common.
#: See the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.CORRELATION_HOVER_FORMAT_FOR_GRAPHS>`__
#: for details.
CORRELATION_HOVER_FORMAT_FOR_GRAPHS: HoverFormat = HoverFormat(jl.MetacellsGraphs.CORRELATION_HOVER_FORMAT_FOR_GRAPHS)

#: The color to give to entities without any value of a category (e.g. with no type or no block). See the Julia
#: `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.EMPTY_CATEGORY_COLOR>`__
#: for details.
EMPTY_CATEGORY_COLOR: str = str(jl.MetacellsGraphs.EMPTY_CATEGORY_COLOR)


def function_default(function: Callable[..., Any], parameter: str) -> Any:
    """
    The default of the ``parameter`` of the Julia function which a ``function`` of this package wraps (e.g.
    ``function_default(fill_n_cells, "log_base")``), as the Julia function declares it. A default which is a Julia
    enum is given as its ``somegraphspy`` enum (e.g. :py:class:`~somegraphspy.LogBase`). A default which depends on
    other parameters (e.g. a title made of the gene's name) is an error. See the Julia
    `documentation <https://tanaylab.github.io/TanayLabUtilities.jl/v0.1.0/documentation.html#TanayLabUtilities.Documentation.function_default>`__
    for details.
    """
    julia_name = function.__name__
    if jl.isdefined(jl.MetacellsGraphs, jl.Symbol(julia_name + "!")):
        julia_name += "!"
    # ``TanayLabUtilities`` is reached through ``MetacellsGraphs``, which uses it.
    value = jl.MetacellsGraphs.DataSources.TanayLabUtilities.function_default(
        jl.getproperty(jl.MetacellsGraphs, jl.Symbol(julia_name)), jl.Symbol(parameter)
    )
    if jl.isa(value, jl.Base.Enum):
        enum_class = getattr(somegraphspy, str(jl.nameof(jl.typeof(value))))
        return enum_class(str(value))
    return value


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
    hover_format: Optional[HoverFormat] = None,
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
        **_given(
            via=_to_julia_array(via),
            empty_value=empty_value,
            entries=_to_julia_array(entries),
            title=title,
            hover_format=hover_format,
        ),
    )


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


# Frame columns.


def fill_column_vector_data(
    sinks: VectorDataSinks,
    frame: pd.DataFrame,
    *,
    column: str,
    title: Optional[str] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Fill the ``sinks`` with the ``column`` of a ``frame``, one value per row. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_column_vector_data!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_column_vector_data_b(
        _to_julia(sinks), _to_julia_frame(frame), column=column, **_given(title=title, hover_format=hover_format)
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
    hover_format: Optional[HoverFormat] = None,
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
        **_given(
            row_entries=_to_julia_array(row_entries),
            column_entries=_to_julia_array(column_entries),
            title=title,
            hover_format=hover_format,
        ),
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


def put_expression_level_configuration(
    sinks: Sinks,
    *,
    gene_fraction_regularization: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
    log_base: Optional[LogBase] = None,
) -> None:
    """
    Show an expression level (of a gene, or of a gene module) on a log scale of the ``log_base`` (2 by default),
    regularized by ``gene_fraction_regularization``, named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_expression_level_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_expression_level_configuration_b(
        _to_julia(sinks),
        **_given(
            gene_fraction_regularization=gene_fraction_regularization,
            title=title,
            show_legend=show_legend,
            log_base=_to_julia(log_base),
        ),
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
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Fill the ``sinks`` with the ``gene`` expression level (linear fraction) per each of the ``entries`` of some ``daf``
    ``axis``, shown by :py:func:`put_expression_level_configuration`, and name the entities after the ``entries``. See
    the Julia
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
            log_base=_to_julia(log_base),
            hover_format=hover_format,
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
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Fill the ``sinks`` with the ``genes`` expression level (linear fraction) per each of the ``entries`` of some ``daf``
    ``axis``, shown by :py:func:`put_expression_level_configuration`, and name the rows and the columns. See the Julia
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
            log_base=_to_julia(log_base),
            hover_format=hover_format,
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
    hover_format: Optional[HoverFormat] = None,
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
            hover_format=hover_format,
        ),
    )


def put_fold_configuration(
    sinks: Sinks, *, max_fold: Optional[float] = None, title: Optional[str] = None, show_legend: Optional[bool] = None
) -> None:
    """
    Show a fold (the log base 2 of a ratio) in the range -``max_fold`` (blue) to 0 (white) to +``max_fold`` (red),
    named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_fold_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_fold_configuration_b(
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
    show_legend: Optional[bool] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Fill the ``sinks`` with the correlation of the ``gene`` per each of the ``entries`` of some ``daf`` ``axis``, shown
    by :py:func:`put_correlation_configuration`, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_gene_correlation!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_gene_correlation_b(
        _to_julia(sinks),
        daf,
        gene=gene,
        **_given(
            axis=axis,
            entries=_to_julia_array(entries),
            title=title,
            show_legend=show_legend,
            hover_format=hover_format,
        ),
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


def get_gene_is_correlated_vector(
    daf: DafReader, *, gene: str, axis: Optional[str] = None, entries: Optional[Entries] = None
) -> np.ndarray:
    """
    Get whether the :py:func:`get_gene_correlation_vector` of the ``gene`` is not zero for each of the ``entries`` of
    some ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_gene_is_correlated_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_gene_is_correlated_vector(
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
    hover_format: Optional[HoverFormat] = None,
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
            axis=axis,
            base_axis=base_axis,
            entries=_to_julia_array(entries),
            title=title,
            show_legend=show_legend,
            hover_format=hover_format,
        ),
    )


def put_correlation_configuration(
    sinks: Sinks,
    values: Union[NumbersVector, NumbersMatrix],
    *,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Show a correlation (or a difference of correlations) on a linear scale, named by the ``title``. As colors, zero is
    white, and the colors follow the range of the ``values``, which always includes zero. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_correlation_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_correlation_configuration_b(
        _to_julia(sinks), _to_julia(values), **_given(title=title, show_legend=show_legend)
    )


def get_correlation_vector(
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
) -> np.ndarray:
    """
    Get the correlation ``property`` of each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_correlation_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_correlation_vector(
            daf,
            axis=axis,
            property=property,
            **_given(entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value),
        )
    )


def fill_correlation(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    title: Optional[str] = None,
    empty_value: Optional[float] = None,
    show_legend: Optional[bool] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Fill the ``sinks`` with the correlation ``property`` of each of the ``entries`` of the ``daf`` ``axis``, shown by
    :py:func:`put_correlation_configuration`, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_correlation!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_correlation_b(
        _to_julia(sinks),
        daf,
        axis=axis,
        property=property,
        **_given(
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            title=title,
            empty_value=empty_value,
            show_legend=show_legend,
            hover_format=hover_format,
        ),
    )


def put_correlation_values(
    sinks: VectorDataSinks,
    value_per_entry: NumbersVector,
    name_per_entry: StringsVector,
    *,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Put the correlation ``value_per_entry`` into the ``sinks``, shown by :py:func:`put_correlation_configuration`, and
    name the entities after the ``name_per_entry``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_correlation_values!>`__
    for details.
    """
    jl.MetacellsGraphs.put_correlation_values_b(
        _to_julia(sinks),
        _to_julia(value_per_entry),
        _to_julia(name_per_entry),
        **_given(title=title, show_legend=show_legend, hover_format=hover_format),
    )


def fill_correlation_between_cells_and_punctuated_metacells(  # pylint: disable=invalid-name
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Fill the ``sinks`` with the correlation of each gene between the cells and their punctuated metacells. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_correlation_between_cells_and_punctuated_metacells!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_correlation_between_cells_and_punctuated_metacells_b(
        _to_julia(sinks),
        daf,
        **_property_fill_given(axis, entries, via, empty_value, title, show_legend, hover_format=hover_format),
    )


def get_correlation_between_cells_and_punctuated_metacells_vector(  # pylint: disable=invalid-name
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
) -> np.ndarray:
    """
    Get the correlation of each gene between the cells and their punctuated metacells. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_correlation_between_cells_and_punctuated_metacells_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_correlation_between_cells_and_punctuated_metacells_vector(
            daf, **_property_get_given(axis, entries, via, empty_value)
        )
    )


def fill_correlation_between_cells_and_projected_punctuated_metacells(  # pylint: disable=invalid-name
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Fill the ``sinks`` with the correlation of each gene between the cells and their projected punctuated metacells. See
    the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_correlation_between_cells_and_projected_punctuated_metacells!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_correlation_between_cells_and_projected_punctuated_metacells_b(
        _to_julia(sinks),
        daf,
        **_property_fill_given(axis, entries, via, empty_value, title, show_legend, hover_format=hover_format),
    )


def get_correlation_between_cells_and_projected_punctuated_metacells_vector(  # pylint: disable=invalid-name
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
) -> np.ndarray:
    """
    Get the correlation of each gene between the cells and their projected punctuated metacells. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_correlation_between_cells_and_projected_punctuated_metacells_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_correlation_between_cells_and_projected_punctuated_metacells_vector(
            daf, **_property_get_given(axis, entries, via, empty_value)
        )
    )


def fill_correlation_between_cells_and_projected_metacells(  # pylint: disable=invalid-name
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Fill the ``sinks`` with the correlation of each gene between the query cells and their projected atlas metacells.
    See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_correlation_between_cells_and_projected_metacells!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_correlation_between_cells_and_projected_metacells_b(
        _to_julia(sinks),
        daf,
        **_property_fill_given(axis, entries, via, empty_value, title, show_legend, hover_format=hover_format),
    )


def get_correlation_between_cells_and_projected_metacells_vector(  # pylint: disable=invalid-name
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
) -> np.ndarray:
    """
    Get the correlation of each gene between the query cells and their projected atlas metacells. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_correlation_between_cells_and_projected_metacells_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_correlation_between_cells_and_projected_metacells_vector(
            daf, **_property_get_given(axis, entries, via, empty_value)
        )
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


def put_coordinates_configuration(sinks: Sinks, *, title: Optional[str] = None) -> None:
    """
    Show the coordinates of an embedding (e.g. a UMAP) without ticks or a grid, named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_coordinates_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_coordinates_configuration_b(_to_julia(sinks), **_given(title=title))


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


# Entry masks.


def get_is_in_categories_vector(
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    categories: Sequence[str],
    via: Optional[Via] = None,
    empty_value: Optional[str] = None,
) -> np.ndarray:
    """
    Get whether the category ``property`` of each entry of the ``daf`` ``axis`` is one of the ``categories``, as a mask
    to pass as the ``entries`` of another data source. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_in_categories_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_in_categories_vector(
            daf,
            axis=axis,
            property=property,
            categories=_to_julia_array(categories),
            **_given(via=_to_julia_array(via), empty_value=empty_value),
        )
    )


def get_is_in_neighborhood_vector(
    daf: DafReader, *, axis: Optional[str] = None, block: str, via: Optional[Via] = None
) -> np.ndarray:
    """
    Get whether each entry of the ``daf`` ``axis`` (by default, the metacells) is in the neighborhood of the ``block``,
    as a mask to pass as the ``entries`` of another data source. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_in_neighborhood_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_in_neighborhood_vector(
            daf, block=block, **_given(axis=axis, via=_to_julia_array(via))
        )
    )


def get_is_in_environment_vector(
    daf: DafReader, *, axis: Optional[str] = None, block: str, via: Optional[Via] = None
) -> np.ndarray:
    """
    Get whether each entry of the ``daf`` ``axis`` (by default, the metacells) is in the environment of the ``block``,
    as a mask to pass as the ``entries`` of another data source. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_in_environment_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_in_environment_vector(daf, block=block, **_given(axis=axis, via=_to_julia_array(via)))
    )


# Counts.


def put_size_configuration(
    sinks: Sinks, *, title: Optional[str] = None, show_legend: Optional[bool] = None, log_base: Optional[LogBase] = None
) -> None:
    """
    Show a size (a count) on a log scale of the ``log_base`` (10 by default), regularized by 1, in the ``YlOrRd`` color
    scale, named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_size_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_size_configuration_b(
        _to_julia(sinks), **_given(title=title, show_legend=show_legend, log_base=_to_julia(log_base))
    )


def get_size_vector(
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[float] = None,
) -> np.ndarray:
    """
    Get the size (count) ``property`` of each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_size_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_size_vector(
            daf,
            axis=axis,
            property=property,
            **_given(entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value),
        )
    )


def fill_size(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    title: Optional[str] = None,
    empty_value: Optional[float] = None,
    show_legend: Optional[bool] = None,
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Fill the ``sinks`` with the size (count) ``property`` of each of the ``entries`` of the ``daf`` ``axis``, shown
    by :py:func:`put_size_configuration`, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_size!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_size_b(
        _to_julia(sinks),
        daf,
        axis=axis,
        property=property,
        **_given(
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            title=title,
            empty_value=empty_value,
            show_legend=show_legend,
            log_base=_to_julia(log_base),
            hover_format=hover_format,
        ),
    )


def put_size_values(
    sinks: VectorDataSinks,
    value_per_entry: NumbersVector,
    name_per_entry: StringsVector,
    *,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
) -> None:
    """
    Put the size (count) ``value_per_entry`` into the ``sinks``, shown by :py:func:`put_size_configuration`, and name
    the entities after the ``name_per_entry``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_size_values!>`__
    for details.
    """
    jl.MetacellsGraphs.put_size_values_b(
        _to_julia(sinks),
        _to_julia(value_per_entry),
        _to_julia(name_per_entry),
        **_given(title=title, show_legend=show_legend, log_base=_to_julia(log_base), hover_format=hover_format),
    )


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
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
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
            log_base=_to_julia(log_base),
            hover_format=hover_format,
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
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
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
            log_base=_to_julia(log_base),
            hover_format=hover_format,
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
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
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
            log_base=_to_julia(log_base),
            hover_format=hover_format,
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
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
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
            axis=axis,
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            title=title,
            show_legend=show_legend,
            log_base=_to_julia(log_base),
            hover_format=hover_format,
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
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
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
            axis=axis,
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            title=title,
            show_legend=show_legend,
            log_base=_to_julia(log_base),
            hover_format=hover_format,
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
    log_base: Optional[LogBase] = None,
    hover_format: Optional[HoverFormat] = None,
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
            axis=axis,
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            title=title,
            show_legend=show_legend,
            log_base=_to_julia(log_base),
            hover_format=hover_format,
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
    colors_axis: Optional[str] = None,
    empty_color: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the block of each of the ``entries`` of the ``daf`` ``axis``, each block in its own color,
    and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_block!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_block_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis,
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            empty_value=empty_value,
            title=title,
            colors_axis=colors_axis,
            empty_color=empty_color,
            show_legend=show_legend,
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


def get_block_colors(
    daf: DafReader, *, empty_value: Optional[str] = None, empty_color: Optional[str] = None
) -> Mapping[str, str]:
    """
    Get the palette of the ``daf`` ``block`` axis, with an additional ``empty_color`` for the ``empty_value`` of entries
    with no block. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_block_colors>`__
    for details.
    """
    colors = jl.MetacellsGraphs.get_block_colors(daf, **_given(empty_value=empty_value, empty_color=empty_color))
    return {str(key): str(jl.getindex(colors, key)) for key in jl.keys(colors)}


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
    via: Optional[Via] = None,
    empty_value: Optional[str] = None,
    title: Optional[str] = None,
    colors_axis: Optional[str] = None,
    empty_color: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the type of each of the ``entries`` of the ``daf`` ``axis``, colored by the palette of the
    ``colors_axis`` (by default, the ``color`` of each type), and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_type!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_type_b(
        _to_julia(sinks),
        daf,
        **_given(
            axis=axis,
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            empty_value=empty_value,
            title=title,
            colors_axis=colors_axis,
            empty_color=empty_color,
            show_legend=show_legend,
        ),
    )


def get_type_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[str] = None,
) -> np.ndarray:
    """
    Get the type of each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_type_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_type_vector(
            daf,
            **_given(axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value),
        )
    )


def get_type_colors(
    daf: DafReader, *, empty_value: Optional[str] = None, empty_color: Optional[str] = None
) -> Mapping[str, str]:
    """
    Get the palette of the ``daf`` ``type`` axis (the ``color`` of each type), with an additional ``empty_color`` for the
    ``empty_value`` of entries with no type. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_type_colors>`__
    for details.
    """
    colors = jl.MetacellsGraphs.get_type_colors(daf, **_given(empty_value=empty_value, empty_color=empty_color))
    return {str(key): str(jl.getindex(colors, key)) for key in jl.keys(colors)}


# Categories.


def put_category_configuration(
    sinks: Sinks, color_per_value: Mapping[str, str], *, title: Optional[str] = None, show_legend: Optional[bool] = None
) -> None:
    """
    Show the values of a category (e.g. types or blocks) of the ``sinks`` using the ``color_per_value`` palette, named
    by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_category_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_category_configuration_b(
        _to_julia(sinks), _to_julia(dict(color_per_value)), **_given(title=title, show_legend=show_legend)
    )


def get_category_colors(
    daf: DafReader, *, axis: str, empty_value: Optional[str] = None, empty_color: Optional[str] = None
) -> Mapping[str, str]:
    """
    Get the palette mapping each entry of the ``daf`` ``axis`` of a category to its color (its ``color`` if the axis
    has one, or a computed one otherwise), with an additional ``empty_color`` for the ``empty_value`` of entries with no
    value. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_category_colors>`__
    for details.
    """
    colors = jl.MetacellsGraphs.get_category_colors(
        daf, axis=axis, **_given(empty_value=empty_value, empty_color=empty_color)
    )
    return {str(key): str(jl.getindex(colors, key)) for key in jl.keys(colors)}


def get_category_vector(
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[str] = None,
) -> np.ndarray:
    """
    Get the category ``property`` (e.g. the type) of each of the ``entries`` of the ``daf`` ``axis``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_category_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_category_vector(
            daf,
            axis=axis,
            property=property,
            **_given(entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value),
        )
    )


def fill_category(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    title: Optional[str] = None,
    empty_value: Optional[str] = None,
    colors_axis: Optional[str] = None,
    empty_color: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the category ``property`` of each of the ``entries`` of the ``daf`` ``axis``, shown by
    :py:func:`put_category_configuration` using the palette of the ``colors_axis`` (by default, the axis named after the
    ``property``), and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_category!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_category_b(
        _to_julia(sinks),
        daf,
        axis=axis,
        property=property,
        **_given(
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            title=title,
            empty_value=empty_value,
            colors_axis=colors_axis,
            empty_color=empty_color,
            show_legend=show_legend,
        ),
    )


def put_category_values(
    sinks: VectorDataSinks,
    value_per_entry: StringsVector,
    name_per_entry: StringsVector,
    color_per_value: Mapping[str, str],
    *,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Put the category ``value_per_entry`` into the ``sinks``, shown by :py:func:`put_category_configuration` using the
    ``color_per_value`` palette, and name the entities after the ``name_per_entry``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_category_values!>`__
    for details.
    """
    jl.MetacellsGraphs.put_category_values_b(
        _to_julia(sinks),
        _to_julia(value_per_entry),
        _to_julia(name_per_entry),
        _to_julia(dict(color_per_value)),
        **_given(title=title, show_legend=show_legend),
    )


# Booleans.


def put_boolean_configuration(sinks: Sinks, *, title: Optional[str] = None, show_legend: Optional[bool] = None) -> None:
    """
    Show Boolean *string* values: ``true`` in black and ``false`` in light grey, named by the ``title``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_boolean_configuration!>`__
    for details.
    """
    jl.MetacellsGraphs.put_boolean_configuration_b(_to_julia(sinks), **_given(title=title, show_legend=show_legend))


def get_boolean_vector(
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get the Boolean ``property`` of each of the ``entries`` of the ``daf`` ``axis``. This can be shown as colors, or
    used as a mask. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_boolean_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_boolean_vector(
            daf,
            axis=axis,
            property=property,
            **_given(entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value),
        )
    )


def fill_boolean(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: str,
    property: str,  # pylint: disable=redefined-builtin
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    title: Optional[str] = None,
    empty_value: Optional[bool] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the Boolean ``property`` of each of the ``entries`` of the ``daf`` ``axis``, shown by
    :py:func:`put_boolean_configuration`, and name the entities after the ``entries``. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_boolean!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_boolean_b(
        _to_julia(sinks),
        daf,
        axis=axis,
        property=property,
        **_given(
            entries=_to_julia_array(entries),
            via=_to_julia_array(via),
            title=title,
            empty_value=empty_value,
            show_legend=show_legend,
        ),
    )


def put_boolean_values(
    sinks: VectorDataSinks,
    value_per_entry: BoolsVector,
    name_per_entry: Optional[StringsVector],
    *,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Put the Boolean ``value_per_entry`` into the ``sinks``, shown by :py:func:`put_boolean_configuration`, and name the
    entities after the ``name_per_entry`` (if given). See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.put_boolean_values!>`__
    for details.
    """
    jl.MetacellsGraphs.put_boolean_values_b(
        _to_julia(sinks),
        _to_julia(value_per_entry),
        _to_julia(name_per_entry),
        **_given(title=title, show_legend=show_legend),
    )


# Gene flags.


def fill_is_lateral(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene is lateral. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_lateral!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_lateral_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_lateral_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene is lateral. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_lateral_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_lateral_vector(daf, **_property_get_given(axis, entries, via, empty_value))
    )


def fill_is_marker(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene is a marker gene. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_marker!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_marker_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_marker_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene is a marker gene. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_marker_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_marker_vector(daf, **_property_get_given(axis, entries, via, empty_value))
    )


def fill_is_transcription_factor(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene is a transcription factor (binds to DNA). See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_transcription_factor!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_transcription_factor_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_transcription_factor_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene is a transcription factor. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_transcription_factor_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_transcription_factor_vector(
            daf, **_property_get_given(axis, entries, via, empty_value)
        )
    )


def fill_is_regulator(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene regulates the expression of other genes. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_regulator!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_regulator_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_regulator_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene regulates the expression of other genes. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_regulator_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_regulator_vector(daf, **_property_get_given(axis, entries, via, empty_value))
    )


def fill_is_forbidden(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene is forbidden. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_forbidden!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_forbidden_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_forbidden_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene is forbidden. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_forbidden_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_forbidden_vector(daf, **_property_get_given(axis, entries, via, empty_value))
    )


def fill_is_skeleton(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene is a skeleton gene. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_skeleton!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_skeleton_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_skeleton_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene is a skeleton gene. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_skeleton_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_skeleton_vector(daf, **_property_get_given(axis, entries, via, empty_value))
    )


def fill_is_correlated_with_skeleton(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene is correlated with some skeleton gene. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_correlated_with_skeleton!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_correlated_with_skeleton_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_correlated_with_skeleton_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene is correlated with some skeleton gene. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_correlated_with_skeleton_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_correlated_with_skeleton_vector(
            daf, **_property_get_given(axis, entries, via, empty_value)
        )
    )


def fill_is_mitochondrial(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene is mitochondrial. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_mitochondrial!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_mitochondrial_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_mitochondrial_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene is mitochondrial. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_mitochondrial_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_mitochondrial_vector(daf, **_property_get_given(axis, entries, via, empty_value))
    )


def fill_is_ribosomal(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene is ribosomal. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_ribosomal!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_ribosomal_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_ribosomal_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene is ribosomal. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_ribosomal_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_ribosomal_vector(daf, **_property_get_given(axis, entries, via, empty_value))
    )


def fill_is_excluded(
    sinks: VectorDataSinks,
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each gene is excluded. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.fill_is_excluded!>`__
    for details.
    """
    jl.MetacellsGraphs.fill_is_excluded_b(
        _to_julia(sinks), daf, **_property_fill_given(axis, entries, via, empty_value, title, show_legend)
    )


def get_is_excluded_vector(
    daf: DafReader,
    *,
    axis: Optional[str] = None,
    entries: Optional[Entries] = None,
    via: Optional[Via] = None,
    empty_value: Optional[bool] = None,
) -> np.ndarray:
    """
    Get whether each gene is excluded. See the Julia
    `documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/data_sources.html#MetacellsGraphs.DataSources.get_is_excluded_vector>`__
    for details.
    """
    return _from_julia_array(
        jl.MetacellsGraphs.get_is_excluded_vector(daf, **_property_get_given(axis, entries, via, empty_value))
    )


# The keywords given to a per-property fill, as passed to Julia.
def _property_fill_given(  # pylint: disable=too-many-positional-arguments
    axis: Optional[str],
    entries: Optional[Entries],
    via: Optional[Via],
    empty_value: Optional[Union[bool, float]],
    title: Optional[str],
    show_legend: Optional[bool],
    *,
    hover_format: Optional[HoverFormat] = None,
) -> Mapping[str, Any]:
    return _given(
        axis=axis,
        entries=_to_julia_array(entries),
        via=_to_julia_array(via),
        empty_value=empty_value,
        title=title,
        show_legend=show_legend,
        hover_format=hover_format,
    )


# The keywords given to a per-property get, as passed to Julia.
def _property_get_given(
    axis: Optional[str], entries: Optional[Entries], via: Optional[Via], empty_value: Optional[Union[bool, float]]
) -> Mapping[str, Any]:
    return _given(axis=axis, entries=_to_julia_array(entries), via=_to_julia_array(via), empty_value=empty_value)


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
