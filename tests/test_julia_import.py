"""
Test that the Julia side is reachable.

Importing this package starts Julia and imports ``MetacellsGraphs``, so a failure here means the environment is not set
up rather than that anything in the package is wrong - which is worth saying separately from any test of what the
package computes.
"""

import json
import os
import tomllib
from typing import Any
from typing import Collection

import metacellsgraphspy as mg
from metacellsgraphspy.julia_import import _GLUE_PROJECT_PATH
from metacellsgraphspy.julia_import import jl


def test_julia_is_reachable() -> None:
    """
    Importing the package starts Julia and imports the Julia package it wraps.
    """
    assert mg.jl_version is not None
    assert mg.jl.MetacellsGraphs is not None


def _glue_project() -> dict[str, Any]:
    """
    The ``Project.toml`` of ``MetacellsGraphsPy``.
    """
    with open(_GLUE_PROJECT_PATH, "rb") as file:
        return tomllib.load(file)


def test_glue_version_is_the_package_version() -> None:
    """
    ``MetacellsGraphsPy`` ships inside this package, so their versions are the same.
    """
    assert _glue_project()["version"] == mg.__version__


def _dependencies_of(package: str) -> Collection[str]:
    """
    The dependencies a Julia package's own ``Project.toml`` lists, if it is loaded in ``Main``.
    """
    if not jl.seval(f"isdefined(Main, :{package})"):
        return ()
    directory = jl.seval(f"pkgdir({package})")
    if directory is None:
        return ()
    with open(os.path.join(str(directory), "Project.toml"), "rb") as file:
        return tomllib.load(file).get("deps", {})


def test_url_packages_are_glue_dependencies() -> None:
    """
    A package comes by URL when it is not registered. It is a dependency of ``MetacellsGraphsPy``, or of one of
    ``MetacellsGraphsPy``'s dependencies.
    """
    with open(os.path.join(os.path.dirname(mg.__file__), "juliapkg.json"), encoding="utf8") as file:
        packages = json.load(file)["packages"]
    glue_packages = _glue_project()["deps"]
    for name, package in packages.items():
        if "url" in package and name != "MetacellsGraphsPy":
            assert name in glue_packages or any(
                name in _dependencies_of(glue_package) for glue_package in glue_packages
            ), f"the URL package {name} of juliapkg.json is not a dependency of MetacellsGraphsPy, nor of its dependencies"
