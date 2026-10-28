"""
Test the calls of a typical session, which ``MetacellsGraphsPy`` makes while it is precompiled.
"""

from metacellsgraphspy.workloads import run_workload


def test_run_workload() -> None:
    """
    The calls succeed, so precompiling ``MetacellsGraphsPy`` succeeds too.
    """
    run_workload()
