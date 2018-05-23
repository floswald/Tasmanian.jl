module Tasmanian

    import Base.show
	using Plots

	# proper build will have to build the library
	const TASlib = "/Applications/TSG/lib/libtasmaniansparsegrid.dylib"

    # includes
    include("TSG.jl")
    include("../examples/examples.jl")

end # module
