"""
Import the Julia environment.

The Julia run-time is obtained by ``dafpy``, which this package depends on, so everything here is taken from
``dafpy.julia_import`` rather than repeated. See it for how Julia is chosen and configured, which is by ``juliacall``'s
own environment variables, plus the ``@default`` value those variables accept here. All of them have to be set before
importing this package, since importing it imports ``dafpy``, which starts Julia.
"""

from dafpy.julia_import import _from_julia_array  # pylint: disable=unused-import
from dafpy.julia_import import _from_julia_frame  # pylint: disable=unused-import
from dafpy.julia_import import _given  # pylint: disable=unused-import
from dafpy.julia_import import _to_julia_array  # pylint: disable=unused-import
from dafpy.julia_import import _to_julia_scalar_or_collection  # pylint: disable=unused-import
from dafpy.julia_import import jl
from dafpy.julia_import import jl_version

__all__ = ["jl", "jl_version"]

# Everything is imported rather than ``using``, so no package's exports leak into Julia's ``Main``. This keeps
# ``Main`` clear for other Python packages that wrap Julia packages and are used in the same session.
for package in ("MetacellsGraphs",):
    jl.seval("import " + package)
