# The Julia side of `metacellsgraphspy`: the Julia packages it uses, and the precompilation of the calls Python makes.
#
# Python always loads `PythonCall` before anything else, so this loads it first too, and is precompiled the same way.
# Its workload runs `metacellsgraphspy.workloads.run_workload()` from Python, so this package caches the code those calls
# compile, and a session does not compile it again.
module MetacellsGraphsPy

using PythonCall
using PrecompileTools

import DafPy
import Metacells
import MetacellsGraphs
import MetacellsPy
import SomeGraphsPy

# The code this caches depends on the Python code of `metacellsgraphspy`, which runs the workload. Julia only knows of its
# own sources, so this tells it to compile again when any of those files changes.
const METACELLSGRAPHSPY_DIRECTORY = joinpath(@__DIR__, "..", "..", "..")
for name in readdir(METACELLSGRAPHSPY_DIRECTORY)
    if endswith(name, ".py")
        include_dependency(joinpath(METACELLSGRAPHSPY_DIRECTORY, name); track_content = true)
    end
end

# Whether the Python that `PythonCall` runs imports the `metacellsgraphspy` this package belongs to. A precompilation
# started from plain Julia gets whatever Python `PythonCall` finds, which may lack `metacellsgraphspy`, or have a version
# of it without this package beside it. A directory named `metacellsgraphspy` in the working directory is also found, as
# a namespace with no origin.
function is_matching_metacellsgraphspy()::Bool
    spec = pyimport("importlib.util").find_spec("metacellsgraphspy")
    if pyis(spec, pybuiltins.None) || pyis(spec.origin, pybuiltins.None)
        return false
    end
    project_path = joinpath(dirname(pyconvert(String, spec.origin)), "julia", "MetacellsGraphsPy", "Project.toml")
    return isfile(project_path) &&
           occursin("uuid = \"$(Base.PkgId(MetacellsGraphsPy).uuid)\"", read(project_path, String))
end

@compile_workload begin
    if !is_matching_metacellsgraphspy()
        error(
            "MetacellsGraphsPy could not be precompiled: the Python that PythonCall runs " *
            "($(PythonCall.C.CTX.exe_path)) does not import the metacellsgraphspy that MetacellsGraphsPy belongs to.\n" *
            "This is harmless: MetacellsGraphsPy will be precompiled the first time you import metacellsgraphspy " *
            "from Python.\n" *
            "To do it now, run: python -c \"import metacellsgraphspy\" " *
            "(with the Python where metacellsgraphspy is installed).",
        )
    end
    # While this is precompiled, neither `metacellsgraphspy` nor the Python packages below it can import anything into
    # Julia's `Main`, so their Julia packages are handed over here instead.
    pysetattr(pyimport("builtins"), "_dafpy_glue", Py(DafPy))
    pysetattr(pyimport("builtins"), "_somegraphspy_glue", Py(SomeGraphsPy))
    pysetattr(pyimport("builtins"), "_metacellspy_glue", Py(MetacellsPy))
    pysetattr(pyimport("builtins"), "_metacellsgraphspy_glue", Py(MetacellsGraphsPy))
    pyexec("from metacellsgraphspy.workloads import run_workload\nrun_workload()", pydict())
end

end  # module
