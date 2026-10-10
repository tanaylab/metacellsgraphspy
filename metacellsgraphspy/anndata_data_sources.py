"""
Data sources of a metacells ``AnnData`` (e.g., read from an ``h5ad`` file).

These are the ``AnnData`` counterparts of some of the ``Daf`` data sources. They are implemented in Python only. Each
``obs`` entry is a metacell, each ``var`` entry is a gene, and ``X`` holds the fraction of the UMIs of each metacell in
each gene. The data they write into the ``somegraphspy`` graph is the same as their ``Daf`` counterparts write.

The ``entries`` are the names or the (1-based) indices of the used metacells, as in the ``Daf`` data sources. As there,
a fill which isn't given its ``entries`` fills the entries its ``sinks`` are already named after, or else all of them.
"""

from typing import Mapping
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp  # type: ignore
from anndata import AnnData  # type: ignore
from somegraphspy import VectorDataSinks
from somegraphspy import get_vector_names_data
from somegraphspy import put_vector_data
from somegraphspy import put_vector_names_data

from .data_sources import EMPTY_TYPE_COLOR
from .data_sources import Entries
from .data_sources import put_genes_expression_configuration
from .data_sources import put_type_configuration

__all__ = [
    "ad_fill_gene_expression",
    "ad_fill_type",
    "ad_get_gene_expression_vector",
    "ad_get_type_colors",
    "ad_get_type_vector",
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
    put_genes_expression_configuration(
        sinks, gene_fraction_regularization=gene_fraction_regularization, title=title, show_legend=show_legend
    )
    put_vector_names_data(sinks, _entry_names(adata, entries))


def ad_get_gene_expression_vector(adata: AnnData, *, gene: str, entries: Optional[Entries] = None) -> np.ndarray:
    """
    Get the vector of the ``gene`` expression level (linear fraction) per each of the ``entries`` of the ``adata``. This
    is the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.get_gene_expression_vector`.
    """
    column = adata.X[:, adata.var_names.get_loc(gene)]
    if sp.issparse(column):
        column = column.toarray()
    return np.asarray(column).ravel()[_entry_indices(adata, entries)]


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
    put_type_configuration(
        sinks,
        ad_get_type_colors(type_colors_csv, empty_type_color=empty_type_color),
        title=title,
        show_legend=show_legend,
    )
    put_vector_names_data(sinks, _entry_names(adata, entries))


def ad_get_type_vector(adata: AnnData, *, type_property: str, entries: Optional[Entries] = None) -> np.ndarray:
    """
    Get the ``type_property`` of the ``obs`` of each of the ``entries`` of the ``adata``. A missing type is returned as
    the empty string. This is the ``AnnData`` counterpart of :py:func:`~metacellsgraphspy.data_sources.get_type_vector`.
    """
    types = adata.obs[type_property].iloc[_entry_indices(adata, entries)]
    return np.array(["" if pd.isna(type_name) else str(type_name) for type_name in types], dtype=str)


def ad_get_type_colors(type_colors_csv: str, *, empty_type_color: Optional[str] = None) -> Mapping[str, str]:
    """
    Get the palette mapping each type to its color from the ``type_colors_csv`` file, with an additional
    ``empty_type_color`` (by default, :py:obj:`~metacellsgraphspy.data_sources.EMPTY_TYPE_COLOR`) for entries with no
    type (the empty string). The file has a header line. Its first column is the type and its second column is the
    color. The names of the columns don't matter. This is the ``AnnData`` counterpart of
    :py:func:`~metacellsgraphspy.data_sources.get_type_colors`.
    """
    frame = pd.read_csv(type_colors_csv, dtype=str, keep_default_na=False)
    color_per_type = dict(zip(frame.iloc[:, 0], frame.iloc[:, 1]))
    color_per_type[""] = EMPTY_TYPE_COLOR if empty_type_color is None else empty_type_color
    return color_per_type


# The entries a fill uses: the given ones, or else the names its sinks already hold. A ``None`` means all of them.
def _sinks_entries(sinks: VectorDataSinks, entries: Optional[Entries]) -> Optional[Entries]:
    if entries is None:
        return get_vector_names_data(sinks)
    return entries


# The 0-based indices of the (1-based or named) entries, or of all the entries.
def _entry_indices(adata: AnnData, entries: Optional[Entries]) -> np.ndarray:
    if entries is None:
        return np.arange(adata.n_obs)
    entries_array = np.asarray(entries)
    if entries_array.dtype.kind in "iu":
        return entries_array - 1
    indices = adata.obs_names.get_indexer(entries_array)
    if np.any(indices < 0):
        raise KeyError(f"unknown entries: {list(entries_array[indices < 0])}")
    return indices


def _entry_names(adata: AnnData, entries: Optional[Entries]) -> np.ndarray:
    return np.asarray(adata.obs_names[_entry_indices(adata, entries)], dtype=str)
