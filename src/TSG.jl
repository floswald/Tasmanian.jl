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

isGlobal(tsg::TasmanianSG)          = convert(Bool,ccall((:tsgIsGlobal,TASlib),Cint,(Ptr{Nothing},),tsg.pGrid))
isSequence(tsg::TasmanianSG)        = convert(Bool,ccall((:tsgIsSequence,TASlib),Cint,(Ptr{Nothing},),tsg.pGrid))
isLocalPolynomial(tsg::TasmanianSG) = convert(Bool,ccall((:tsgIsLocalPolynomial,TASlib),Cint,(Ptr{Nothing},),tsg.pGrid))
isWavelet(tsg::TasmanianSG)         = convert(Bool,ccall((:tsgIsWavelet,TASlib),Cint,(Ptr{Nothing},),tsg.pGrid))

getNumDimensions(tsg::TasmanianSG) = convert(Int,ccall((:tsgGetNumDimensions,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))
getDims(tsg::TasmanianSG)          = tsg.dims
getNumLoaded(tsg::TasmanianSG)     = convert(Int,ccall((:tsgGetNumLoaded,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))
getNumPoints(tsg::TasmanianSG)     = convert(Int,ccall((:tsgGetNumPoints,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))
getNumOutputs(tsg::TasmanianSG)    = convert(Int,ccall((:tsgGetNumOutputs,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))
getNout(tsg::TasmanianSG)          = tsg.nout
getNumNeeded(tsg::TasmanianSG)     = convert(Int,ccall((:tsgGetNumNeeded,TASlib),Int32,(Ptr{Nothing},),tsg.pGrid))

# make polynomial grid 
function makeLocalPolynomialGrid!(TSG::TasmanianSG; iOrder::Int=1, sRule::AbstractString="localp", ilevelLimits=Int[])

    n = length(ilevelLimits)
    levelLimits = C_NULL
    if n > 0
        if n != iDimension
            throw(ArgumentError("invalid number of levellimites. levelLimits = $levelLimits needs to have iDimension elements"))
        else
            levelLimits = ilevelLimits
        end
    end

    ccall((:tsgMakeLocalPolynomialGrid,TASlib),Nothing,(Ptr{Nothing},Cint,Cint,Cint,Cint,Cstring,Ptr{Nothing}),TSG.pGrid, TSG.dims, TSG.nout, TSG.depth, iOrder, sRule, levelLimits)

end

function setSurplusRefinement!(TSG::TasmanianSG, tol::Float64; iOutput::Int=-1, sCriteria::AbstractString="", ilevelLimits=Int[])

    if getNumLoaded(TSG) == 0
        throw(MethodError("cannot call setSurplusRefinement for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    if tol < 0
        throw(ArgumentError("tol needs to be a positive number"))
    end
    n = length(ilevelLimits)
    levelLimits = C_NULL
    if n > 0
        if n != iDimension
            throw(ArgumentError("invalid number of levellimites. levelLimits = $levelLimits needs to have iDimension elements"))
        else
            levelLimits = ilevelLimits
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
            ccall((:tsgSetLocalSurplusRefinement,TASlib),Nothing,(Ptr{Nothing},Cdouble,Cstring,Cint,Ptr{Nothing}),TSG.pGrid,tol,sCriteria,iOutput,levelLimits)
        end
    end
end

function getPoints(TSG::TasmanianSG)
    iNumDims = getNumDimensions(TSG)
    iNumPoints = getNumPoints(TSG)
    if iNumPoints == 0
        zeros(2)
    else
        out = zeros(Float64,iNumPoints*iNumDims)
        ccall((:tsgGetPointsStatic,TASlib),Nothing,(Ptr{Nothing},Ptr{Cdouble}),TSG.pGrid,out)
    end
    return reshape(out,iNumDims,iNumPoints)'
end

function getNeededPoints(TSG::TasmanianSG)
    iNumDims = getNumDimensions(TSG)
    iNumPoints = getNumNeeded(TSG)
    if iNumPoints == 0
        return zeros(2)
    else
        out = zeros(Float64,iNumPoints*iNumDims)
        ccall((:tsgGetNeededPointsStatic,TASlib),Nothing,(Ptr{Nothing},Ptr{Cdouble}),TSG.pGrid,out)
    end
    return reshape(out,iNumDims,iNumPoints)'
end



# load needed points
function loadNeededPoints!(tsg::TasmanianSG,vals::Array{Float64})
    n = size(vals)
    if length(n) > 1
        if n[1] != getNumNeeded(tsg)
            if getNumNeeded(tsg) == 0
                if n[1] != getNumLoaded(tsg)
                    throw(ArgumentError("vals does not have $(getNumLoaded(tsg)) rows"))
                end
            else
                throw(ArgumentError("vals does not have $(getNumNeeded(tsg)) rows"))
            end
        elseif n[2] != getNumOutputs(tsg)
            throw(ArgumentError("vals does not have $(getNumOutputs(tsg)) columns"))
        end
        # vals = vals[:]
    else
        if n[1] != getNumNeeded(tsg)
            if getNumNeeded(tsg) == 0
                if n[1] != getNumLoaded(tsg)
                    throw(ArgumentError("vals does not have $(getNumLoaded(tsg)) rows"))
                end
            else
                throw(ArgumentError("vals does not have $(getNumNeeded(tsg)) rows"))
            end
        end
    end

    ccall((:tsgLoadNeededPoints,TASlib),Nothing,(Ptr{Nothing},Ptr{Cdouble}),tsg.pGrid,vals)
end

function evaluateBatch(tsg::TasmanianSG,vals::Matrix{Float64})
    if getNumLoaded(tsg) == 0
        throw(MethodError("cannot call evaluateBatch for a grid before any points are loaded, i.e., call loadNeededPoints first!"))
    end
    n = size(vals)
    nx = n[1]
    nd = n[2]
    if nd != getDims(tsg)
        throw(ArgumentError("vals does not have $(getNumOutputs(tsg)) columns"))
    end
    v = reshape(vals',nx*nd,1)
    out = zeros(Float64,nx*tsg.nout)

    ccall((:tsgEvaluateBatch,TASlib),Nothing,(Ptr{Nothing},Ptr{Cdouble},Cint,Ptr{Cdouble}),tsg.pGrid,v,nx,out)
    return reshape(out,nx,tsg.nout)

end

function setDomainTransform!(tsg::TasmanianSG,doms::Matrix{Float64})
	n,m = size(doms)
	if n != getNumDimensions(tsg)
		throw(ArgumentError("doms must have `getNumDimensions(tsg)`=$(getNumDimensions(tsg)) rows"))
	end
	if m != 2
		throw(ArgumentError("doms must have 2 columns"))
	end
	pA = doms[:,1]
	pB = doms[:,2]
	ccall((:tsgSetDomainTransform,TASlib),Nothing,(Ptr{Nothing},Ptr{Cdouble},Ptr{Cdouble}),tsg.pGrid,pA,pB)
end

    # def evaluateBatch(self, llfX):
    # '''
    # evaluates the intepolant at the points of interest and returns
    # the result

    # this should be called after the grid has been created and after
    # values have been loaded

    # llfX: a 2-D numpy.ndarray
    #       with second dimension equal to iDimensions
    #       each row in the array is a single requested point

    # output: a 2-D numpy.ndarray
    #         with dimensions llfX.shape[0] X iOutputs
    #         each row corresponds to the value of the interpolant
    #         for one row of llfX

    # '''
    # if (self.getNumLoaded() == 0):
    #     raise TasmanianInputError("evaluateBatch", "ERROR: cannot call evaluateBatch for a grid before any points are loaded, i.e., call loadNeededPoints first!")
    # if (len(llfX.shape) != 2):
    #     raise TasmanianInputError("llfX", "ERROR: llfX should be a 2-D numpy.ndarray instread it has dimension {0:1d}".format(len(llfX.shape)))
    # iNumX = llfX.shape[0]
    # if (iNumX == 0):
    #     return np.empty([0, self.getNumOutputs()], np.float64)
    # iNumDim = llfX.shape[1]
    # if (iNumDim != self.getNumDimensions()):
    #     raise TasmanianInputError("llfX", "ERROR: llfX.shape[1] should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumDim))
    # iNumOutputs = self.getNumOutputs()
    # aY = np.empty([iNumX, iNumOutputs], np.float64)
    # # np.ctypeslib.as_ctypes(llfX.reshape([iNumX*iNumDim,])) messes up, the first 4 entries randomly get set to machine eps (10^-310) and 0
    # lfX = llfX.reshape([iNumX*iNumDim,])
    # self.pLibTSG.tsgEvaluateBatch(self.pGrid, np.ctypeslib.as_ctypes(lfX), iNumX, np.ctypeslib.as_ctypes(aY.reshape([iNumX*iNumOutputs,])))
    # return aY

# def loadNeededPoints(self, llfVals):
#     '''
#     loads the values of the target function at the needed points
#     if there are no needed points, this reset the currently loaded
#     values

#     llfVals: a 2-D numpy.ndarray
#              with dimensions getNumNeeded() X iOutputs
#              each row corresponds to the values of the outputs at
#              the corresponding needed point. The order and leading
#              dimension must match the points obtained form
#              getNeededPoints()

#     '''
#     if (len(llfVals.shape) != 2):
#         raise TasmanianInputError("llfVals", "ERROR: llfVals should be a 2-D numpy.ndarray, instead it has {0:1d} dimensions".format(len(llfVals.shape)))
#     if (llfVals.shape[0] != self.getNumNeeded()):
#         if (self.getNumNeeded() == 0):
#             if (llfVals.shape[0] != self.getNumLoaded()):
#                 raise TasmanianInputError("llfVals", "ERROR: leading dimension of llfVals is {0:1d} but the number of current points is {1:1d}".format(llfVals.shape[0], self.getNumNeeded()))
#         else:
#             raise TasmanianInputError("llfVals", "ERROR: leading dimension of llfVals is {0:1d} but the number of needed points is {1:1d}".format(llfVals.shape[0], self.getNumNeeded()))
#     if (llfVals.shape[1] != self.getNumOutputs()):
#         raise TasmanianInputError("llfVals", "ERROR: second dimension of llfVals is {0:1d} but the number of outputs is set to {1:1d}".format(llfVals.shape[1], self.getNumOutputs()))
#     iNumPoints = llfVals.shape[0]
#     iNumDims = llfVals.shape[1]
#     self.pLibTSG.tsgLoadNeededPoints(self.pGrid, np.ctypeslib.as_ctypes(llfVals.reshape([iNumPoints * iNumDims])))

# getters
# -------
# def getAlpha(self):
#     '''
#     returns the value of fAlpha in the call to makeGlobalGrid
#     if makeGlobalGrid has not been called, returns 0.0

#     '''
#     return self.pLibTSG.tsgGetAlpha(self.pGrid)

# def getBeta(self):
#     '''
#     returns the value of fBeta in the call to makeGlobalGrid
#     if makeGlobalGrid has not been called, returns 0.0

#     '''
#     return self.pLibTSG.tsgGetBeta(self.pGrid)

# def getOrder(self):
#     '''
#     returns the value of iOrder in the call to
#     makeLocalPolynomialGrid or makeWaveletGrid

#     if makeLocalPolynomialGrid and makeWaveletGrid
#     have not been called, returns -1

#     '''
#     return self.pLibTSG.tsgGetOrder(self.pGrid)


# def getNumOutputs(self):
#     '''
#     returns the value of iOutputs in the make***Grid command
#     if no grid has been made, it returns 0

#     '''
#     return self.pLibTSG.tsgGetNumOutputs(self.pGrid)

# def getRule(self):
#     '''
#     returns the value of sRule in the make***Grid command
#     if makeWaveletGrid is used, returns "wavelet"
#     if no grid has been made, it returns "unknown"

#     '''
#     sRule = self.pLibTSG.tsgGetRule(self.pGrid)
#     if (sys.version_info.major == 3):
#         sRule = str(sRule, encoding='utf8')
#     return sRule

# def getCustomRuleDescription(self):
#     '''
#     returns the description provided in the custom rule file
#     if not using a custom grid, returns ""

#     '''
#     if ("custom-tabulated" in self.getRule()):
#         sRule = self.pLibTSG.tsgGetCustomRuleDescription(self.pGrid)
#         if (sys.version_info.major == 3):
#             sRule = str(sRule, encoding='utf8')
#         return sRule
#     else:
#         return ""

# def getNumNeeded(self):
#     '''
#     returns the number of points needed to form the interpolant or
#     form the next interpolant following a refinement

#     '''
#     return self.pLibTSG.tsgGetNumNeeded(self.pGrid)


# def getLoadedPoints(self):
#     '''
#     returns the points loaded in the existing interpolant

#     output: a 2-D numpy.ndarray of size getNumLoaded() X iDimension
#         reach row correspoinds to one point
#         if (getNumLoaded() == 0): returns numpy.empty([0,0])

#     '''
#     iNumDims = self.getNumDimensions()
#     iNumPoints = self.getNumLoaded()
#     if (iNumPoints == 0):
#         return np.empty([0, 0], np.float64)
#     aPoints = np.empty([iNumPoints * iNumDims], np.float64)
#     self.pLibTSG.tsgGetLoadedPointsStatic(self.pGrid, np.ctypeslib.as_ctypes(aPoints))
#     return aPoints.reshape([iNumPoints, iNumDims])

# def getNeededPoints(self):
#     '''
#     returns the points needed to form the interpolant or the next
#     level of refinement following a set***Refinement() call

#     output: 2-D numpy.ndarray of size getNumNeeded() X iDimension
#         reach row correspoinds to one point
#         if (getNumNeeded() == 0): returns numpy.empty([0,0])

#     '''
#     iNumDims = self.getNumDimensions()
#     iNumPoints = self.getNumNeeded()
#     if (iNumPoints == 0):
#         return np.empty([0, 0], np.float64)
#     aPoints = np.empty([iNumPoints * iNumDims], np.float64)
#     self.pLibTSG.tsgGetNeededPointsStatic(self.pGrid, np.ctypeslib.as_ctypes(aPoints))
#     return aPoints.reshape([iNumPoints, iNumDims])


# def getPoints(self):
#     '''
#     if points have been loaded, gives the same as getLoadedPoints()
#     otherwise, returns the same as getNeededPoints()

#     '''
#     iNumDims = self.getNumDimensions()
#     iNumPoints = self.getNumPoints()
#     if (iNumPoints == 0):
#         return np.empty([0, 0], np.float64)
#     aPoints = np.empty([iNumPoints * iNumDims], np.float64)
#     self.pLibTSG.tsgGetPointsStatic(self.pGrid, np.ctypeslib.as_ctypes(aPoints))
#     return aPoints.reshape([iNumPoints, iNumDims])

# def getQuadratureWeights(self):
#     '''
#     returns the quadrature weights associated with
#     the points in getPoints()

#     output: a 1-D numpy.ndarray of length getNumPoints()
#             the order of the weights matches
#             the order in getPoints()

#     '''
#     iNumPoints = self.getNumPoints()
#     if (iNumPoints == 0):
#         return np.empty([0], np.float64)
#     aWeights = np.empty([iNumPoints], np.float64)
#     self.pLibTSG.tsgGetQuadratureWeightsStatic(self.pGrid, np.ctypeslib.as_ctypes(aWeights))
#     return aWeights

# def getInterpolationWeights(self, lfX):
#     '''
#     returns the interpolation weights associated with the points
#     in getPoints()

#     lfX: a 1-D numpy.ndarray with length iDimensions
#          the entries indicate the points for evaluating the weights

#     output: a 1-D numpy.ndarray of length getNumPoints()
#         the order of the weights matches the order in getPoints()

#     '''
#     iNumX = len(lfX)
#     if (iNumX != self.getNumDimensions()):
#         raise TasmanianInputError("lfX", "ERROR: len(lfX) should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumX))
#     iNumPoints = self.getNumPoints()
#     if (iNumPoints == 0):
#         return np.empty([0], np.float64)
#     aWeights = np.empty([iNumPoints], np.float64)
#     self.pLibTSG.tsgGetInterpolationWeightsStatic(self.pGrid, np.ctypeslib.as_ctypes(lfX), np.ctypeslib.as_ctypes(aWeights))
#     return aWeights

# def getInterpolationWeightsBatch(self, llfX):
#     '''
#     returns the interpolation weights associated with the points
#     in getPoints()

#     finds multiple weights with a single library call
#     uses OpenMP if enabled in libtasmaniansparsegrids.so

#     llfX: a 2-D numpy.ndarray with second dimension iDimensions
#           each row in the array is a single requested point

#     output: a 2-D numpy.ndarray
#             with dimensions llfX.shape[0] X getNumPoints()
#             each row corresponds to the weight for one row of llfX

#     '''
#     if (len(llfX.shape) != 2):
#         raise TasmanianInputError("llfX", "ERROR: llfX should be a 2-D numpy.ndarray instread it has dimension {0:1d}".format(len(llfX.shape)))
#     iNumX = llfX.shape[0]
#     if (iNumX == 0):
#         return np.empty([0, self.getNumPoints()], np.float64)
#     iNumDim = llfX.shape[1]
#     if (iNumDim != self.getNumDimensions()):
#         raise TasmanianInputError("llfX", "ERROR: llfX.shape[1] should equal {0:1d} instead it equals {1:1d}".format(self.getNumDimensions(), iNumDim))
#     iNumPoints = self.getNumPoints()
#     if (iNumPoints == 0):
#         return np.empty([0, 0], np.float64)
#     aWeights = np.empty([iNumX, iNumPoints], np.float64)
#     self.pLibTSG.tsgBatchGetInterpolationWeightsStatic(self.pGrid,
#         np.ctypeslib.as_ctypes(llfX.reshape([iNumX * iNumDim])), iNumX,
#         np.ctypeslib.as_ctypes(aWeights.reshape([iNumX * iNumPoints])))
#     return aWeights
