LocalRules = ["localp" "localp-zero" "semi-localp" "localp-boundary"]


mutable struct TasmanianSG
    pGrid   :: Ptr{Nothing}
    version :: VersionNumber
    dims    :: Int
    nout    :: Int
    depth   :: Int

	function TasmanianSG(dims::Int,nout::Int,depth::Int)
		this = new()
		output_ptr = ccall(
	        (:tsgConstructTasmanianSparseGrid,TASlib), # name of C function and library
	        Ptr{TasmanianSG},              # output type
	        ()                        # tuple of input types
		)
	    if output_ptr == C_NULL # Could not allocate memory
	        throw(OutOfMemoryError())
	    else
	    	this.pGrid = output_ptr
	    end
        this.version = VersionNumber(unsafe_string(ccall((:tsgGetVersion,TASlib),Cstring,())))
        this.dims    = dims
        this.nout    = nout
        this.depth   = depth
	    return this
	end
end

function show(io::IO,TSG::TasmanianSG)
    ccall((:tsgPrintStats,TASlib),Nothing,(Ptr{Nothing},),TSG.pGrid)
end
# function del(TSG::TasmanianSG)
#   ccall((:tsgDestructTasmanianSparseGrid,TASlib),Nothing,(Ptr{Nothing},),TSG.pGrid)
# end

"""
    isGlobal(tsg::TasmanianSG)

returns True if using a global grid
"""
isGlobal(tsg::TasmanianSG)          = convert(Bool,ccall((:tsgIsGlobal,TASlib),Cint,(Ptr{Nothing},),tsg.pGrid))

"""
    isSequence(tsg::TasmanianSG)

returns true if using a sequence grid
"""
isSequence(tsg::TasmanianSG)        = convert(Bool,ccall((:tsgIsSequence,TASlib),Cint,(Ptr{Nothing},),tsg.pGrid))

"""
    isLocalPolynomial(tsg::TasmanianSG)

returns true if using a local polynomial grid
"""
isLocalPolynomial(tsg::TasmanianSG) = convert(Bool,ccall((:tsgIsLocalPolynomial,TASlib),Cint,(Ptr{Nothing},),tsg.pGrid))

"""
    isWavelet(tsg::TasmanianSG)
returns true if using a local wavelet grid
"""
isWavelet(tsg::TasmanianSG)         = convert(Bool,ccall((:tsgIsWavelet,TASlib),Cint,(Ptr{Nothing},),tsg.pGrid))

"""
    isFourier(tsg::TasmanianSG)

returns true if using a Fourier grid
"""
isFourier(tsg::TasmanianSG)         = convert(Bool,ccall((:tsgIsFourier,TASlib),Cint,(Ptr{Nothing},),tsg.pGrid))

"""
    getNumDimensions(tsg::TasmanianSG)

returns the value of iDimension in the make***Grid command
if no grid has been made, it returns 0
"""
getNumDimensions(tsg::TasmanianSG) = convert(Int,ccall((:tsgGetNumDimensions,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))

getDims(tsg::TasmanianSG)          = tsg.dims

"""
    getNumLoaded(tsg::TasmanianSG)

returns the number of points loaded in the existing interpolant
"""
getNumLoaded(tsg::TasmanianSG)     = convert(Int,ccall((:tsgGetNumLoaded,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))

"""
    getNumNeeded(tsg::TasmanianSG)

returns the number of points needed to form the interpolant or
        form the next interpolant following a refinement
"""
getNumNeeded(tsg::TasmanianSG)     = convert(Int,ccall((:tsgGetNumNeeded,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))

"""
    getNumPoints(tsg::TasmanianSG)

if points have been loaded, returns the same as getNumLoaded()
otherwise, returns the same as getNumNeeded()
"""
getNumPoints(tsg::TasmanianSG)     = convert(Int,ccall((:tsgGetNumPoints,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))

"""
    getNumOutputs(tsg::TasmanianSG)

returns the value of iOutputs in the make***Grid command
if no grid has been made, it returns 0
"""
getNumOutputs(tsg::TasmanianSG)    = convert(Int,ccall((:tsgGetNumOutputs,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))

getNout(tsg::TasmanianSG)          = tsg.nout

"""
    getOrder(tsg::TasmanianSG)

returns the value of iOrder in the call to
makeLocalPolynomialGrid or makeWaveletGrid
if makeLocalPolynomialGrid and makeWaveletGrid
have not been called, returns -1
"""
getOrder(tsg::TasmanianSG)          = convert(Int,ccall((:tsgGetOrder,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))

# make polynomial grid 
"""
    makeLocalPolynomialGrid!(TSG::TasmanianSG; iOrder::Int=1, sRule::AbstractString="localp", lilevelLimits=Int[])

creates a new sparse grid using a local polynomial rule
discards any existing grid held by TSG

iDimension: int (positive)
      the number of inputs (default = 1)

iOutputs: int (non-negative)
      the number of outputs

iDepth: int (non-negative)
        controls the density of the grid, i.e.,
        the number of levels to use

iOrder: int (must be -1 or bigger)
        -1 indicates largest possible order
         1 means linear, 2 means quadratic, etc.
         0 means piece-wise constant, it has different hierarchy
           then the other orders, most notably the 1D rule
           triples the number of points per level (as opposed
           to double for the other cases)

sRule: string (defines the 1-D rule that induces the grid)
      'localp' (default) 'localp-zero'  'semi-localp'  'localp-boundary'

"""
function makeLocalPolynomialGrid!(TSG::TasmanianSG; iOrder::Int=1, sRule::AbstractString="localp", lilevelLimits=Int[])
    if iOrder < -1
        throw(ArgumentError("order should be a non-negative integer"))
    end
    if !(sRule in LocalRules)
        throw(ArgumentError("invalid local polynomial rule, see TasmanianSG.LocalRules for list of accepted sequence rules"))
    end
    
    n = length(lilevelLimits)
    levelLimits = C_NULL
    if n > 0
        if n != iDimension
            throw(ArgumentError("invalid number of levellimites. lilevelLimits = $lilevelLimits needs to have $iDimension elements"))
        else
            levelLimits = lilevelLimits
        end
    end

    ccall((:tsgMakeLocalPolynomialGrid,TASlib),Nothing,(Ptr{Nothing},Cint,Cint,Cint,Cint,Cstring,Ptr{Nothing}),TSG.pGrid, TSG.dims, TSG.nout, TSG.depth, iOrder, sRule, levelLimits)
end

"""
    setSurplusRefinement!(TSG::TasmanianSG, tol::Float64; iOutput::Int=-1, sCriteria::AbstractString="", lilevelLimits=Int[], llfScaleCorrection = Float64[])

using hierarchical surplusses as an error indicator, the surplus
refinement adds points to the grid to improve accuracy

when using sequence grids: this algorithm corresponds to the
                           greedy Knapsack problem

when using local polynomial or wavelet grids, this call
                         corresponds to local spatial refinement

fTolerance: float (non-negative)
            the relative error tolerance, i.e.,
            we refine only for points associated with surplus
            that exceeds the tolerance

iOutput: int (indicates the output to use)
         selects which output to use for refinement
         sequence and local polynomial grids accept -1 to
         indicate all outputs

sCriteria: hierarchical and direction refinement strategy
           'classic'  'parents'   'direction'   'fds'   'stable'
          applicable only for Local Polynomial and Wavelet grids

llfScaleCorrection: 2-D numpy.ndarray of non-negative numbers
                    Instead of comparing the normalized surpluses to
                    the tolerance, the scaled surplus will be used.
                    The correction allows to manually guide the
                    refinement process.
                    The surplus of the j-th output of the i-th point
                    will be scaled by llfScaleCorrection[i][j].
                    llfScaleCorrection.shape[0] must be equal to
                    getNumLoaded()
                    If empty, the scale is assumed 1.0
                    llfScaleCorrection.shape[1] must be equal to the
                    number of outputs used in the process, which is
                    equal to getNumOutputs() for iOutput == -1,
                    or 1 if iOutput > -1.
"""
function setSurplusRefinement!(TSG::TasmanianSG, tol::Float64; iOutput::Int=-1, sCriteria::AbstractString="", lilevelLimits=Int[], llfScaleCorrection = Float64[])

    if getNumLoaded(TSG) == 0
        throw(MethodError("cannot call setSurplusRefinement for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    if tol < 0
        throw(ArgumentError("tol needs to be a positive number"))
    end
    n = length(lilevelLimits)
    levelLimits = C_NULL
    if n > 0
        if n != iDimension
            throw(ArgumentError("invalid number of levellimites. lilevelLimits = $lilevelLimits needs to have $iDimension elements"))
        else
            levelLimits = lilevelLimits
        end
    end
    if length(sCriteria)==0
        if !isSequence(TSG)
            throw(MethodError("sCriteria must be specified"))
        else
            ccall((:tsgSetGlobalSurplusRefinement,TASlib),Nothing,(Ptr{Nothing},Cdouble,Cint,Ptr{Nothing}),TSG.pGrid,tol,iOutput,levelLimits)
        end
    else
        if isSequence(TSG)
            throw(MethodError("sCriteria not used for Sequence Grids"))
        else
            ccall((:tsgSetLocalSurplusRefinement,TASlib),Nothing,(Ptr{Nothing},Cdouble,Cstring,Cint,Ptr{Nothing},Ptr{Nothing}),TSG.pGrid,tol,sCriteria,iOutput,levelLimits, llfScaleCorrection)
        end
    end
end

"""
    getPoints(TSG::TasmanianSG)

if points have been loaded, gives the same as getLoadedPoints()
otherwise, returns the same as getNeededPoints()
"""
function getPoints(TSG::TasmanianSG)
    iNumDims = getNumDimensions(TSG)
    iNumPoints = getNumPoints(TSG)
    if iNumPoints == 0
        zeros(2)
    else
        out = zeros(Float64, iNumDims, iNumPoints)
        ccall((:tsgGetPointsStatic,TASlib), Nothing, (Ptr{Nothing}, Ptr{Cdouble}), TSG.pGrid, out)
    end
    return out
end

"""
    getLoadedPoints(TSG::TasmanianSG)

returns the points loaded in the existing interpolant

output: 2-D array of size iDimension X getNumNeeded()
    each column  corresponds to one point
    if (getNumNeeded() == 0): returns numpy.empty([0,0])
"""
function getLoadedPoints(TSG::TasmanianSG)
    iNumDims = getNumDimensions(TSG)
    iNumPoints = getNumLoaded(TSG)
    if iNumPoints == 0
        return zeros(2)
    else
        out = zeros(Float64, iNumDims, iNumPoints)
        ccall((:tsgGetNeededPointsStatic, TASlib), Nothing, (Ptr{Nothing}, Ptr{Cdouble}), TSG.pGrid, out)
    end
    return out
end

"""
    getNeededPoints(TSG::TasmanianSG)

returns the points needed to form the interpolant or the next
level of refinement following a set***Refinement() call

output: 2-D array of size iDimension X getNumNeeded()
    each column  corresponds to one point
    if (getNumNeeded() == 0): returns an array iDmension X 0
"""
function getNeededPoints(TSG::TasmanianSG)
    iNumDims = getNumDimensions(TSG)
    iNumPoints = getNumNeeded(TSG)
    if iNumPoints == 0
        return Matrix{Float64}(undef, iNumDims, 0)
    else
        out = zeros(Float64, iNumDims, iNumPoints)
        ccall((:tsgGetNeededPointsStatic, TASlib), Nothing, (Ptr{Nothing}, Ptr{Cdouble}), TSG.pGrid, out)
    end
    return out
end

# load needed points
"""
    loadNeededPoints!(tsg::TasmanianSG, vals::Array{Float64})

loads the values of the target function at the needed points
if there are no needed points, this reset the currently loaded
values

vals: an array with dimensions iOutputs X getNumNeeded() 
      each column corresponds to the values of the outputs at
      the corresponding needed point. The order and leading
      dimension must match the points obtained from
      getNeededPoints()
"""
function loadNeededPoints!(tsg::TasmanianSG, vals::Array{Float64})
    numOutputs = getNumOutputs(tsg)
    nd = ndims(vals)
    if nd == 1
        n = length(vals)
        if mod(n, numOutputs) == 0
            n1 = numOutputs
            n2 = length(vals)/numOutputs
        else
            throw(ArgumentError("vals is a vector but its length isn't a multiple of the numOutputs"))
        end
    elseif nd == 2
        n1, n2 = size(vals)
    else
        throw(ArgumentError("vals must be a vector or a matrix"))
    end
    
    if n1 != numOutputs
        throw(ArgumentError("leading dimension of vals is $n1 but the number of outputs is set to $(getNumOutputs(tsg))"))
    end
    if getNumNeeded(tsg) == 0
        if n2 != getNumLoaded(tsg)
            throw(ArgumentError("the second dimension of vals is $n2 but the number of current points is $(getNumLoaded(tsg))"))
        end
    elseif n2 != getNumNeeded(tsg)
        throw(ArgumentError("the second dimension of vals is $n2 but the number of needed points is $(getNumNeeded(tsg))"))
    end
    ccall((:tsgLoadNeededPoints,TASlib),Nothing,(Ptr{Nothing},Ptr{Cdouble}),tsg.pGrid,vals)
end

"""
     evaluateBatch(tsg::TasmanianSG, vals::Matrix{Float64})

evaluates the intepolant at the points of interest and returns
the result

this should be called after the grid has been created and after
values have been loaded

vals: a 1-D or 2-D array
      with first dimension equal to iDimensions
      each column in the array is a single requested point

output: a 1-D or 2-D matrix
        with dimensions iOutputs X size(vals, 2) 
        each columns corresponds to the value of the interpolant
        for one columns of vals
"""
function evaluateBatch(tsg::TasmanianSG, vals::VecOrMat{Float64})
    iNumOutputs = getNumOutputs(tsg)
    iNumX = size(vals, 2)
    if iNumX > 1
        aY = Matrix{Float64}(undef, iNumOutputs, iNumX)
    else
        aY = Vector{Float64}(undef, iNumOutputs)
    end 
    evaluateBatch!(aY, tsg, vals)
    return aY
end

"""
     evaluateBatch!(aY::VecOrMat{Float64}, tsg::TasmanianSG, vals::VecOrMat{Float64})

evaluates the intepolant at the points of interest and set the result in `aY`

this should be called after the grid has been created and after
values have been loaded

aY: a 1-D or 2-D matrix
        with dimensions iOutputs X size(vals, 2) 
        each columns corresponds to the value of the interpolant
        for one columns of vals

vals: a 1-D or 2-D array
      with first dimension equal to iDimensions
      each column in the array is a single requested point

"""
function evaluateBatch!(aY::VecOrMat{Float64}, tsg::TasmanianSG, vals::VecOrMat{Float64})
    if getNumLoaded(tsg) == 0
        throw(ArgumentError("cannot call evaluateBatch for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    n1 = size(vals)
    if ndims(vals) > 2
        throw(ArgumentError("vals should be a 1-D or 2-D array instead it has $(len(n)) dimensions"))
    end
    @assert ndims(aY) == ndims(vals) "aY and vals must have the same number of columns"
    iNumDim = n1[1]
    iNumX = (length(n1) == 2) ? n1[2] : 1
    n2 = size(aY)
    iNumDimOut = n2[1]
    iNumXOut = (length(n2) == 2) ? n2[2] : 1
    @assert iNumX == iNumXOut "aY and vals must have the same number of columns"
    iNumX == 0 && return
    if (iNumDim != getNumDimensions(tsg))
        throw(ArgumentError("ERROR: size(vals, 1) should equal $(getNumDimensions(tsg)) instead it equals $iNumDim"))
    end
    if (iNumDimOut != getNumOutputs(tsg))
        throw(ArgumentError("ERROR: size(aY, 1) should equal $(getNumOutputs(tsg)) instead it equals $iNumDimOut"))
    end
    ccall((:tsgEvaluateBatch,TASlib),Nothing,(Ptr{Nothing},Ptr{Cdouble},Cint,Ptr{Cdouble}), tsg.pGrid, vals, iNumX, aY)
    return aY
end

"""
    setDomainTransform!(tsg::TasmanianSG, Transform::Matrix{Float64})

sets the lower and upper bound for each dimension

Note: gauss-aguerre and gauss-hermite rules are defined on
      unbounded domain, in which case  this sets the
      shift and scale parameters, consult the manual

Transform: a 2-D array of size iDimension X 2
           transform specifies the lower and upper bound
           of the domain in each direction.

           For gauss-laguerre and gauss-hermite grids, the
           transform gives the a and b parameters of the
           weights
            exp(-b (x - a))
            exp(-b (x - a)^2)
"""
function setDomainTransform!(tsg::TasmanianSG, Transform::Matrix{Float64})
    n = size(Transform)
    if length(n) != 2
        throw(ArgumentError("Transform should be a 2-D array"))
    end
    if n[1] != getNumDimensions(tsg)
        throw(ArgumentError("the first dimension of Transform is $(n[1]) and it should match iDimension: $(getNumDimensions(tsg))"))
    end
    if n[2] != 2
        throw(ArgumentError("the second dimension of Transform is $(n[2]) and it should be 2"))
    end
    pA = Transform
    pB = view(Transform, :, 2)
    ccall((:tsgSetDomainTransform,TASlib),Nothing,(Ptr{Nothing},Ptr{Cdouble},Ptr{Cdouble}),tsg.pGrid,pA,pB)
end

"""
    copyGrid(tsg::TasmanianSG, iOutputsBegin = 0, iOutputsEnd = -1)

accepts an instance of TasmanianSparseGrid class and creates
a hard copy of the class and all included data
original class is not modified

pGrid: instance of TasmanianSparseGrid class
    the source for the copy

iOutputsBegin: integer indicating the first output to copy

iOutputsEnd: integer one bigger than the last output to copy
             if set to -1, all outputs from iOutputsBegin to
             the end will be copied

Examples:

grid.copyGrid(other, 0, -1) # copy all outputs (default)
grid.copyGrid(other, 0, other.getNumOutputs()) # also copy all
grid.copyGrid(other, 0, 3) # copy outputs 0, 1, and 2
grid.copyGrid(other, 1, 4) # copy outputs 1, 2, and 3
"""
function copyGrid(tsg::TasmanianSG, iOutputsBegin = 0, iOutputsEnd = -1)
    newgrid = TasmanianSG(tsg.dims, tsg.nout, tsg.depth)
    ccall((:tsgCopySubGrid,TASlib), Nothing, (Ptr{Nothing}, Ptr{Nothing}, Cint, Cint), newgrid.pGrid, tsg.pGrid, iOutputsBegin, iOutputsEnd)
    return newgrid
end


    
