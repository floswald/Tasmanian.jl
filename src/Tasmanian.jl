module Tasmanian

	if VERSION >= v"0.7.0-DEV.3382"
    	import Libdl
	end

	# Load in `deps.jl`, complaining if it does not exist
	const depsjl_path = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
	if !isfile(depsjl_path)
	    error("Tasmanian.jl not installed properly, run Pkg.build(\"Tasmanian\"), restart Julia and try again")
	end
	include(depsjl_path)

	# Module initialization function
	function __init__()
	    # Always check your dependencies from `deps.jl`
	    check_deps()
	end

    import Base.show
	using Plots

    # includes
    include("TSG.jl")
    include("../examples/examples.jl")

end # module
