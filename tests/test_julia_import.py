"""
Test that the Julia side is reachable.

Importing this package starts Julia and imports ``MetacellsGraphs``, so a failure here means the environment is not set
up rather than that anything in the package is wrong - which is worth saying separately from any test of what the
package computes.
"""

import metacellsgraphspy as mg


def test_julia_is_reachable() -> None:
    """
    Importing the package starts Julia and imports the Julia package it wraps.
    """
    assert mg.jl_version is not None
    assert mg.jl.MetacellsGraphs is not None
