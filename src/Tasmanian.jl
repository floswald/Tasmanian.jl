module Tasmanian

import Base.show
using Plots
using Random
using Tasmanian_jll

const global TASlib = libtasmaniansparsegrid

# includes
include("TSG.jl")
include("../examples/examples.jl")

export TasmanianSG, LocalRules, copyGrid, evaluateBatch, evaluateBatch!, getDims,
getLoadedPoints, getNeededPoints, getNout, getNumDimensions, getNumLoaded,
getNumNeeded, getNumOutputs, getNumPoints, getOrder, getPoints, isFourier,
isGlobal, isLocalPolynomial, isSequence, isWavelet, loadNeededPoints!,
makeLocalPolynomialGrid!, setDomainTransform!, setSurplusRefinement!

end # module
