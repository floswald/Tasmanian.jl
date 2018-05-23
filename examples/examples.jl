    function run()
        iDim = 2
        iOut = 1
        iDepth = 5
        tsg = Tasmanian.TasmanianSG(iDim,iOut,iDepth)
        println(tsg.version)
        println(tsg.pGrid)
        which_basis = 1 #1= linear basis functions -> Check the manual for other options
        Tasmanian.makeLocalPolynomialGrid!(tsg,iOrder=which_basis,sRule="localp")
        return tsg
    end

    unitbox(x) = x*2 - 1

    function ex1()
        iDim = 2
        iOut = 1
        iDepth = 5
        tsg = Tasmanian.TasmanianSG(iDim,iOut,iDepth)
        which_basis = 1 #1= linear basis functions -> Check the manual for other options
        Tasmanian.makeLocalPolynomialGrid!(tsg,iOrder=which_basis,sRule="localp")

        # sparse grid points from that object
        spPoints = getPoints(tsg)

        # measure perf at N randomly chosen points
        tfun(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y)
        srand(1)
        N = 1000
        randPnts = unitbox.(rand(N,2))

        # truth
        truth = [tfun(randPnts[i,1],randPnts[i,2]) for i in 1:N]

        # values on sparse grid
        spVals = [tfun(spPoints[i,1],spPoints[i,2]) for i in 1:size(spPoints,1)]

        # load points needed for such values
        loadNeededPoints!(tsg,spVals)

        # evaluate interpolation
        res = evaluateBatch(tsg,randPnts)

        # Plots
        p1=scatter(randPnts[:,1],randPnts[:,2],m=(:black,2),title="test points")
        p2=scatter(spPoints[:,1],spPoints[:,2],m=(:red,2),title="sparse grid")
        pts=plot(p1,p2,legend=false)

        p4 = surface(randPnts[:,1],randPnts[:,2],truth.-res,title="error")
        p3 = surface(randPnts[:,1],randPnts[:,2],res,title="prediction")
        p = plot(p1,p2,p3,p4,layout=(2,2),legend=false)

        # # compute error 
        return (p,maximum(abs,res .- truth))
    end

    function ex2()
        dim =  2       
        outs = 1       
        iDepth = 2
        tol = 1e-5    
        K = 7  # max refinement steps
        tsg = Tasmanian.TasmanianSG(dim,outs,iDepth)
        which_basis = 1 #1= linear basis functions -> Check the manual for other options
        Tasmanian.makeLocalPolynomialGrid!(tsg,iOrder=which_basis,sRule="localp")

        # sparse grid points from that object
        spPoints = getPoints(tsg)

        # test fun 
        tfun(x,y) = exp(-x^2) * cos(y)

        srand(2)
        N = 1000
        randPnts = unitbox.(rand(N,2))
        # truth
        truth = [tfun(randPnts[i,1],randPnts[i,2]) for i in 1:N]

        # values on sparse grid
        spVals = [tfun(spPoints[i,1],spPoints[i,2]) for i in 1:size(spPoints,1)]
        # load points needed for such values
        loadNeededPoints!(tsg,spVals)

        # evaluate interpolation
        res = evaluateBatch(tsg,randPnts)

        numpoints = size(spPoints,1)

        info("error on initial grid:    $(round(maximum(abs,res .- truth),5)), with $numpoints points")

        # refinefment loop
        anim = @animate for k in 1:K
            setSurplusRefinement!(tsg,tol,sCriteria="classic")
            spPoints = getNeededPoints(tsg)   # additional set of points required after refinement
            spVals = [tfun(spPoints[i,1],spPoints[i,2]) for i in 1:size(spPoints,1)]
            # load points needed for such values
            loadNeededPoints!(tsg,spVals)
            numpoints =+ size(spPoints,1)

            # evaluate interpolation
            res = evaluateBatch(tsg,randPnts)
            info("refinement level $k error: $(round(maximum(abs,res .- truth),5)), with $numpoints points")
        	plot(scatter(spPoints[:,1],spPoints[:,2],title="level $k grid: $numpoints points",m=(:black,1,:+)),
        		 scatter3d(randPnts[:,1],randPnts[:,2],res[:],title="level $k error: $(round(maximum(abs,res .- truth),5))",m=(:red,1)),
        		 layout=(1,2),leg=false
        		 )
        end
        gif(anim,joinpath(dirname(@__FILE__),"ex2.gif"),fps=1)
    end