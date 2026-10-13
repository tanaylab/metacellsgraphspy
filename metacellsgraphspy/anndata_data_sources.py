"""
Data sources of a metacells ``AnnData`` (e.g., read from an ``h5ad`` file).

These are the ``AnnData`` counterparts of some of the ``Daf`` data sources. They are implemented in Python only. Each
``obs`` entry is a metacell, each ``var`` entry is a gene, and ``X`` holds the fraction of the UMIs of each metacell in
each gene. The data they write into the ``somegraphspy`` graph is the same as their ``Daf`` counterparts write.

Each data source reads a fixed column of the Python metacells layout: ``x`` and ``y`` (the UMAP), and ``total_umis``
of the ``obs``, and ``lateral_gene``, ``excluded_gene`` and ``marker_gene`` of the ``var``. The type is the exception,
since its column is named differently in many data sets, so it is given as the ``type_property``.

The ``entries`` are the names, the (1-based) indices, or a Boolean mask of the used metacells (or, for a property of
the genes, of the used genes), as in the ``Daf`` data sources. As there, a fill which isn't given its ``entries`` fills
the entries its ``sinks`` are already named after, or else all of them.
"""

from typing import Mapping
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp  # type: ignore
from anndata import AnnData  # type: ignore
from somegraphspy import LogBase
from somegraphspy import VectorDataSinks
from somegraphspy import get_vector_names_data
from somegraphspy import put_vector_data
from somegraphspy import put_vector_names_data

from .data_sources import EMPTY_CATEGORY_COLOR
from .data_sources import Entries
from .data_sources import put_boolean_values
from .data_sources import put_category_configuration
from .data_sources import put_coordinates_configuration
from .data_sources import put_expression_level_configuration
from .data_sources import put_size_values
from .data_sources import put_umap_data

__all__ = [
    "ad_fill_gene_expression",
    "ad_fill_is_excluded",
    "ad_fill_is_lateral",
    "ad_fill_is_marker",
    "ad_fill_total_UMIs",
    "ad_fill_type",
    "ad_fill_umap",
    "ad_get_gene_expression_vector",
    "ad_get_is_excluded_vector",
    "ad_get_is_lateral_vector",
    "ad_get_is_marker_vector",
    "ad_get_total_UMIs_vector",
    "ad_get_type_colors",
    "ad_get_type_vector",
    "ad_get_umap_vector",
]


def ad_fill_gene_expression(
    sinks: VectorDataSinks,
    adata: AnnData,
    *,
    gene: str,
    entries: Optional[Entries] = None,
    gene_fraction_regularization: Optional[float] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with the ``gene`` expression level (linear fraction) per each of the ``entries`` of the
    ``adata``, shown in log base 2, and name the entities after the ``entries``. The default ``title`` is
    ``<gene> fraction``. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.data_sources.fill_gene_expression`.
    """
    if title is None:
        title = f"{gene} fraction"
    entries = _sinks_entries(sinks, entries)
    put_vector_data(sinks, ad_get_gene_expression_vector(adata, gene=gene, entries=entries), title=title)
    put_expression_level_configuration(
        sinks, gene_fraction_regularization=gene_fraction_regularization, title=title, show_legend=show_legend
    )
    put_vector_names_data(sinks, _entry_names(adata.obs_names, entries))


def ad_get_gene_expression_vector(adata: AnnData, *, gene: str, entries: Optional[Entries] = None) -> np.ndarray:
    """
    Get the vector of the ``gene`` expression level (linear fraction) per each of the ``entries`` of the ``adata``. This
    is the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.get_gene_expression_vector`.
    """
    column = adata.X[:, adata.var_names.get_loc(gene)]
    if sp.issparse(column):
        column = column.toarray()
    return np.asarray(column).ravel()[_entry_indices(adata.obs_names, entries)]


def ad_fill_type(
    sinks: VectorDataSinks,
    adata: AnnData,
    *,
    type_property: str,
    type_colors_csv: str,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
    empty_type_color: Optional[str] = None,
) -> None:
    """
    Fill the ``sinks`` with the ``type_property`` of each of the ``entries`` of the ``adata``, colored by the
    ``type_colors_csv`` file, and name the entities after the ``entries``. The default ``title`` is ``type``. See
    :py:func:`ad_get_type_vector` and :py:func:`ad_get_type_colors` for details. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.data_sources.fill_type`.
    """
    if title is None:
        title = "type"
    entries = _sinks_entries(sinks, entries)
    put_vector_data(sinks, ad_get_type_vector(adata, type_property=type_property, entries=entries), title=title)
    put_category_configuration(
        sinks,
        ad_get_type_colors(type_colors_csv, empty_type_color=empty_type_color),
        title=title,
        show_legend=show_legend,
    )
    put_vector_names_data(sinks, _entry_names(adata.obs_names, entries))


def ad_get_type_vector(adata: AnnData, *, type_property: str, entries: Optional[Entries] = None) -> np.ndarray:
    """
    Get the ``type_property`` of the ``obs`` of each of the ``entries`` of the ``adata``. A missing type is returned as
    the empty string. This is the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.get_type_vector`.
    """
    types = adata.obs[type_property].iloc[_entry_indices(adata.obs_names, entries)]
    return np.array(["" if pd.isna(type_name) else str(type_name) for type_name in types], dtype=str)


def ad_get_type_colors(type_colors_csv: str, *, empty_type_color: Optional[str] = None) -> Mapping[str, str]:
    """
    Get the palette mapping each type to its color from the ``type_colors_csv`` file, with an additional
    ``empty_type_color`` (by default, :py:obj:`~metacellsgraphspy.data_sources.EMPTY_CATEGORY_COLOR`) for entries with no
    type (the empty string). The file has a header line. Its first column is the type and its second column is the
    color. The names of the columns don't matter. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.data_sources.get_type_colors`.
    """
    frame = pd.read_csv(type_colors_csv, dtype=str, keep_default_na=False)
    color_per_type = dict(zip(frame.iloc[:, 0], frame.iloc[:, 1]))
    color_per_type[""] = EMPTY_CATEGORY_COLOR if empty_type_color is None else empty_type_color
    return color_per_type


def ad_fill_umap(
    sinks: VectorDataSinks,
    adata: AnnData,
    *,
    coordinate: str,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
) -> None:
    """
    Fill the ``sinks`` with the UMAP ``coordinate`` (``x`` or ``y``) of each of the ``entries`` of the ``adata``, and
    name the entities after the ``entries``. The default ``title`` is ``UMAP X`` (or ``UMAP Y``). This is the
    ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.fill_umap`.
    """
    if title is None:
        title = f"UMAP {coordinate.upper()}"
    entries = _sinks_entries(sinks, entries)
    put_umap_data(sinks, ad_get_umap_vector(adata, coordinate=coordinate, entries=entries), title=title)
    put_coordinates_configuration(sinks, title=title)
    put_vector_names_data(sinks, _entry_names(adata.obs_names, entries))


def ad_get_umap_vector(adata: AnnData, *, coordinate: str, entries: Optional[Entries] = None) -> np.ndarray:
    """
    Get the UMAP ``coordinate`` (the ``x`` or ``y`` of the ``obs``) of each of the ``entries`` of the ``adata``. This is
    the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.get_umap_vector`.
    """
    return np.asarray(adata.obs[coordinate], dtype="float32")[_entry_indices(adata.obs_names, entries)]


def ad_fill_total_UMIs(  # pylint: disable=invalid-name
    sinks: VectorDataSinks,
    adata: AnnData,
    *,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
    log_base: Optional[LogBase] = None,
) -> None:
    """
    Fill the ``sinks`` with the total UMIs of each of the ``entries`` of the ``adata``, shown as a size, and name the
    entities after the ``entries``. The default ``title`` is ``total UMIs``. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.data_sources.fill_total_UMIs`.
    """
    if title is None:
        title = "total UMIs"
    entries = _sinks_entries(sinks, entries)
    put_size_values(
        sinks,
        ad_get_total_UMIs_vector(adata, entries=entries),
        _entry_names(adata.obs_names, entries),
        title=title,
        show_legend=show_legend,
        log_base=log_base,
    )


def ad_get_total_UMIs_vector(  # pylint: disable=invalid-name
    adata: AnnData, *, entries: Optional[Entries] = None
) -> np.ndarray:
    """
    Get the total UMIs (the ``total_umis`` of the ``obs``) of each of the ``entries`` of the ``adata``. This is the
    ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.get_total_UMIs_vector`.
    """
    return np.asarray(adata.obs["total_umis"])[_entry_indices(adata.obs_names, entries)]


def ad_fill_is_lateral(
    sinks: VectorDataSinks,
    adata: AnnData,
    *,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each of the gene ``entries`` of the ``adata`` is lateral, and name the entities after
    the ``entries``. The default ``title`` is ``lateral``. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.data_sources.fill_is_lateral`.
    """
    if title is None:
        title = "lateral"
    _fill_var_flag(sinks, adata, column="lateral_gene", entries=entries, title=title, show_legend=show_legend)


def ad_get_is_lateral_vector(adata: AnnData, *, entries: Optional[Entries] = None) -> np.ndarray:
    """
    Get whether each of the gene ``entries`` of the ``adata`` is lateral (the ``lateral_gene`` of the ``var``). This is
    the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.get_is_lateral_vector`.
    """
    return _var_flag_vector(adata, column="lateral_gene", entries=entries)


def ad_fill_is_excluded(
    sinks: VectorDataSinks,
    adata: AnnData,
    *,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each of the gene ``entries`` of the ``adata`` is excluded, and name the entities
    after the ``entries``. The default ``title`` is ``excluded``. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.data_sources.fill_is_excluded`.
    """
    if title is None:
        title = "excluded"
    _fill_var_flag(sinks, adata, column="excluded_gene", entries=entries, title=title, show_legend=show_legend)


def ad_get_is_excluded_vector(adata: AnnData, *, entries: Optional[Entries] = None) -> np.ndarray:
    """
    Get whether each of the gene ``entries`` of the ``adata`` is excluded (the ``excluded_gene`` of the ``var``). This
    is the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.get_is_excluded_vector`.
    """
    return _var_flag_vector(adata, column="excluded_gene", entries=entries)


def ad_fill_is_marker(
    sinks: VectorDataSinks,
    adata: AnnData,
    *,
    entries: Optional[Entries] = None,
    title: Optional[str] = None,
    show_legend: Optional[bool] = None,
) -> None:
    """
    Fill the ``sinks`` with whether each of the gene ``entries`` of the ``adata`` is a marker gene, and name the
    entities after the ``entries``. The default ``title`` is ``marker``. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.data_sources.fill_is_marker`.
    """
    if title is None:
        title = "marker"
    _fill_var_flag(sinks, adata, column="marker_gene", entries=entries, title=title, show_legend=show_legend)


def ad_get_is_marker_vector(adata: AnnData, *, entries: Optional[Entries] = None) -> np.ndarray:
    """
    Get whether each of the gene ``entries`` of the ``adata`` is a marker gene (the ``marker_gene`` of the ``var``).
    This is the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.get_is_marker_vector`.
    """
    return _var_flag_vector(adata, column="marker_gene", entries=entries)


# Fill the sinks with a Boolean ``column`` of the ``var``, named after the gene ``entries``.
def _fill_var_flag(
    sinks: VectorDataSinks,
    adata: AnnData,
    *,
    column: str,
    entries: Optional[Entries],
    title: str,
    show_legend: Optional[bool],
) -> None:
    entries = _sinks_entries(sinks, entries)
    put_boolean_values(
        sinks,
        _var_flag_vector(adata, column=column, entries=entries),
        _entry_names(adata.var_names, entries),
        title=title,
        show_legend=show_legend,
    )


# A Boolean ``column`` of the ``var`` for each of the gene ``entries``.
def _var_flag_vector(adata: AnnData, *, column: str, entries: Optional[Entries]) -> np.ndarray:
    return np.asarray(adata.var[column], dtype=bool)[_entry_indices(adata.var_names, entries)]


# The entries a fill uses: the given ones, or else the names its sinks already hold. A ``None`` means all of them.
def _sinks_entries(sinks: VectorDataSinks, entries: Optional[Entries]) -> Optional[Entries]:
    if entries is None:
        return get_vector_names_data(sinks)
    return entries


# The 0-based indices of the (1-based, named or masked) entries among the ``names`` of an axis, or of all of them.
def _entry_indices(names: pd.Index, entries: Optional[Entries]) -> np.ndarray:
    if entries is None:
        return np.arange(len(names))
    entries_array = np.asarray(entries)
    if entries_array.dtype.kind == "b":
        if len(entries_array) != len(names):
            raise ValueError(f"mask length: {len(entries_array)} is different from the number of entries: {len(names)}")
        return np.flatnonzero(entries_array)
    if entries_array.dtype.kind in "iu":
        return entries_array - 1
    indices = names.get_indexer(pd.Index(entries_array))
    if np.any(indices < 0):
        raise KeyError(f"unknown entries: {list(entries_array[indices < 0])}")
    return indices


# The names of the entries among the ``names`` of an axis.
def _entry_names(names: pd.Index, entries: Optional[Entries]) -> np.ndarray:
    return np.asarray(names[_entry_indices(names, entries)], dtype=str)
