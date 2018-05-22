module Tasmanian

	# using Plots

	# proper build will have to build the library
	const TASlib = "/Applications/TSG/lib/libtasmaniansparsegrid.dylib"

	mutable struct TasmanianSG
		pGrid :: Ptr{Void}
		version :: VersionNumber

		function TasmanianSG()
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
		    return this
		end
	end

# 	function gethostname()
#     hostname = Vector{UInt8}(128)
#     ccall((:gethostname, "libc"), Int32,
#           (Ptr{UInt8}, Csize_t),
#           hostname, sizeof(hostname))
#     hostname[end] = 0; # ensure null-termination
#     return unsafe_string(pointer(hostname))
# end



	# methods
	# -------

	# version

	# isOpenMP

	# read/write from disk

	# make globalgrid


	# make sequence grid

	# make polynomial grid 
	function makeLocalPolynomialGrid!(TSG::TasmanianSG,iDimension::Int,iOutputs::Int,iDepth::Int;iOrder=1,sRule="localp",levelLimits=Int[])

		n = length(levelLimits)
		if n > 0
			if n != iDimension
				throw(ArgumentError("invalid number of levellimites. levelLimits = $levelLimits is invalid. must be equal to iDimension"))
			end
		end

		ccall((:tsgMakeLocalPolynomialGrid,TASlib),Void,(Ptr{Void},Cint,Cint,Cint,Cint,Ptr{UInt8},Ptr{Cint}),TSG.pGrid, iDimension, iOutputs, iDepth, iOrder, sRule, levelLimits)

	end

	# make Wavelet grid 

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

    getNumDimensions(tsg::TasmanianSG) = ccall((:tsgGetNumDimensions,TASlib),Int32,(Ptr{Void},),tsg.pGrid)

    # def getNumDimensions(self):
    #     '''
    #     returns the value of iDimension in the make***Grid command
    #     if no grid has been made, it returns 0

    #     '''
    #     return self.pLibTSG.tsgGetNumDimensions(self.pGrid)

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

    getNumLoaded(tsg::TasmanianSG) = ccall((:tsgGetNumLoaded,TASlib),Int32,(Ptr{Void},),tsg.pGrid)

    # def getNumLoaded(self):
    #     '''
    #     returns the number of points loaded in the existing interpolant

    #     '''
    #     return self.pLibTSG.tsgGetNumLoaded(self.pGrid)

    # def getNumNeeded(self):
    #     '''
    #     returns the number of points needed to form the interpolant or
    #     form the next interpolant following a refinement

    #     '''
    #     return self.pLibTSG.tsgGetNumNeeded(self.pGrid)

    getNumPoints(tsg::TasmanianSG) = ccall((:tsgGetNumPoints,TASlib),Int32,(Ptr{Void},),tsg.pGrid)

    # def getNumPoints(self):
    #     '''
    #     if points have been loaded, returns the same as getNumLoaded()
    #     otherwise, returns the same as getNumNeeded()

    #     '''
    #     return self.pLibTSG.tsgGetNumPoints(self.pGrid)

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

    function getPoints(TSG::TasmanianSG)
    	iNumDims = getNumDimensions(TSG)
    	iNumPoints = getNumPoints(TSG)
    	if iNumPoints == 0
    		zeros(2)
    	else
    		out = zeros(Float64,iNumPoints,iNumDims)
    		ccall((:tsgGetPointsStatic,TASlib),Void,(Ptr{Void},Ptr{Cdouble}),TSG.pGrid,out)
    	end
    	return out
    end

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



    function run()
    	tsg = Tasmanian.TasmanianSG()
    	println(tsg.version)
    	makeLocalPolynomialGrid!(tsg,2,1,1)
    	println("num dims = $(getNumDimensions(tsg))")
    	println("num points = $(getNumPoints(tsg))")
    	println("num loaded = $(getNumLoaded(tsg))")
    	println("points:")
    	getPoints(tsg)
    end






end # module
