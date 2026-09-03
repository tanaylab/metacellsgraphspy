"""
Draw the figures worth looking at of a metacells repository.

This is a thin wrapper for the ``MetacellsGraphs.jl`` Julia package, which sits between ``somegraphspy``, which draws
from arrays and knows nothing of metacells, and ``metacellspy``, which computes and draws nothing. A figure reads a
repository, does the arithmetic the picture needs, and returns something a notebook cell can display.

The API works "just the same" as the Julia one, so the documentation mostly links to the relevant entry in the Julia
`documentation <https://tanaylab.github.io/MetacellsGraphs.jl/v0.1.0/index.html>`__.
"""

__author__ = "Oren Ben-Kiki"
__email__ = "oren@ben-kiki.org"
__version__ = "0.1.0"

# pylint: disable=wildcard-import,unused-wildcard-import

from .heatmap_graphs import *
from .julia_import import *
from .scatter_graphs import *
