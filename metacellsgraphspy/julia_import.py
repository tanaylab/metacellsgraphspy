"""
Import the Julia environment.

The Julia run-time is obtained by ``dafpy``, which this package depends on, so everything here is taken from
``dafpy.julia_import`` rather than repeated. See it for how Julia is chosen and configured, which is by ``juliacall``'s
own environment variables, plus the ``@default`` value those variables accept here. All of them have to be set before
importing this package, since importing it imports ``dafpy``, which starts Julia.
"""

import builtins
import os
import tomllib
from typing import Any

from dafpy.julia_import import DEFAULT  # pylint: disable=unused-import
from dafpy.julia_import import DefaultValue  # pylint: disable=unused-import
from dafpy.julia_import import _from_julia_array  # pylint: disable=unused-import
from dafpy.julia_import import _from_julia_frame  # pylint: disable=unused-import
from dafpy.julia_import import _given  # pylint: disable=unused-import
from dafpy.julia_import import _to_julia_array  # pylint: disable=unused-import
from dafpy.julia_import import _to_julia_frame  # pylint: disable=unused-import
from dafpy.julia_import import _to_julia_scalar_or_collection  # pylint: disable=unused-import
from dafpy.julia_import import jl as _dafpy_jl
from dafpy.julia_import import jl_version

__all__ = ["jl", "jl_version"]

# The Julia package ``MetacellsGraphsPy``, in ``julia/MetacellsGraphsPy`` beside this file, caches the code our calls
# compile. Its ``Project.toml`` lists the Julia packages we use. Like ``metacellspy``, this package needs no Julia code of
# its own.
_GLUE_PROJECT_PATH = os.path.join(os.path.dirname(__file__), "julia", "MetacellsGraphsPy", "Project.toml")

# While ``MetacellsGraphsPy`` itself is precompiled, Julia may not import anything into ``Main``, so
# ``MetacellsGraphsPy`` hands itself over instead, and it binds the same names.
_METACELLSGRAPHSPY_GLUE: Any = getattr(builtins, "_metacellsgraphspy_glue", None)

#: The interface to the Julia run-time. This is Julia's ``Main``, as in ``dafpy``, except while ``MetacellsGraphsPy``
#: itself is precompiled.
jl = _dafpy_jl if _METACELLSGRAPHSPY_GLUE is None else _METACELLSGRAPHSPY_GLUE

if _METACELLSGRAPHSPY_GLUE is None:
    with open(_GLUE_PROJECT_PATH, "rb") as glue_project_file:
        _GLUE_PACKAGES = tomllib.load(glue_project_file)["deps"]
    # Everything is imported rather than ``using``, so no package's exports leak into Julia's ``Main``. This keeps
    # ``Main`` clear for other Python packages that wrap Julia packages and are used in the same session.
    for package in ["MetacellsGraphsPy"] + sorted(_GLUE_PACKAGES):
        if package != "PrecompileTools":
            jl.seval("import " + package)
