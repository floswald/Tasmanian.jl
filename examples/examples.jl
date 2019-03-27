    function run()
        iDim = 2
        iOut = 1
        iDepth = 5
        tsg = Tasmanian.TasmanianSG(iDim,iOut,iDepth)
        which_basis = 1 #1= linear basis functions -> Check the manual for other options
        Tasmanian.makeLocalPolynomialGrid!(tsg,iOrder=which_basis,sRule="localp")
        return tsg
    end

    negbox(x) = x*2 - 1
    neg2unit(x) = (x+1)/2

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
        Random.seed!(1)
        N = 1000
        randPnts = negbox.(rand(N,2))

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
        return 0
    end

    function ex2(;save=false)
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

        Random.seed!(2)
        N = 1000
        randPnts = negbox.(rand(N,2))
        # truth
        truth = [tfun(randPnts[i,1],randPnts[i,2]) for i in 1:N]

        # values on sparse grid
        spVals = [tfun(spPoints[i,1],spPoints[i,2]) for i in 1:size(spPoints,1)]
        # load points needed for such values
        loadNeededPoints!(tsg,spVals)

        # evaluate interpolation
        res = evaluateBatch(tsg,randPnts)

        numpoints = size(spPoints,1)

        @info("error on initial grid:    $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")

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
	        pred = Tasmanian.evaluateBatch(tsg,Array(spPoints))  # prediction on spGrid
            @info("refinement level $k error: $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")

            # plot
            zerone = (-1.1,1.1)
        	plot(scatter(spPoints[:,1],spPoints[:,2],title="level $k grid:\n $numpoints points",m=(:black,1,:+),aspect_ratio=:equal,xlims=zerone,ylims=zerone),
        		 scatter3d(randPnts[:,1],randPnts[:,2],res[:],title="max error: $(round(maximum(abs,res .- truth),digits = 5))",m=(:red,1),xlims=zerone,ylims=zerone,zlims=(0,1)),
        		 scatter3d(spPoints[:,1],spPoints[:,2],pred[:],title="grid prediction",m=(:red,1,0.2),zlims=(0,1),xlims=zerone,ylims=zerone,zgrid=:black),
        		 layout=(1,3),leg=false
        		 )
        end
        if save 
            gif(anim,joinpath(dirname(@__FILE__),"ex2.gif"),fps=1)
        end
        return 0
    end


    function ex3(;save=false)
    	tfun(x,y) = 1.0 / (abs(0.5 - x^4 - y^4) + 0.1)
        dim =  2       
        outs = 1       
        iDepth = 2
        tol = 1e-5    
        K = 9  # max refinement steps
        tsg = Tasmanian.TasmanianSG(dim,outs,iDepth)

        which_basis = 1 #1= linear basis functions -> Check the manual for other options
        Tasmanian.makeLocalPolynomialGrid!(tsg,iOrder=which_basis,sRule="localp")

        # domain is [0,1] here
        setDomainTransform!(tsg,[0 1.0;0 1])

        # sparse grid points from that object
        spPoints = Tasmanian.getPoints(tsg)
        Random.seed!(2)
        N = 1000
        randPnts = rand(N,2)
        # truth
        truth = [tfun(randPnts[i,1],randPnts[i,2]) for i in 1:N]

        # values on sparse grid
        spVals = [tfun(spPoints[i,1],spPoints[i,2]) for i in 1:size(spPoints,1)]
        # load points needed for such values
        loadNeededPoints!(tsg,spVals)

        # evaluate interpolation
        res = evaluateBatch(tsg,randPnts)
        # res = evaluateBatch(tsg,randPnts)

        numpoints = size(spPoints,1)

        @info("error on initial grid:    $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")

        # refinefment loop
        anim = @animate for k in 1:K
            Tasmanian.setSurplusRefinement!(tsg,tol,sCriteria="classic")
	        spPoints = Tasmanian.getNeededPoints(tsg)
            spVals = [tfun(spPoints[i,1],spPoints[i,2]) for i in 1:size(spPoints,1)]
            # load points needed for such values
            Tasmanian.loadNeededPoints!(tsg,spVals)
            numpoints =+ size(spPoints,1)

            # evaluate interpolation
	        res = Tasmanian.evaluateBatch(tsg,randPnts)
	        pred = Tasmanian.evaluateBatch(tsg,Array(spPoints))  # prediction on spGrid
	        # res = evaluateBatch(tsg,randPnts)
	        zerone = (-0.1,1.1)
            @info("refinement level $k error: $(round(maximum(abs,res .- truth),digits = 5)), with $numpoints points")
        	plot(scatter(spPoints[:,1],spPoints[:,2],title="level $k grid: $numpoints points",m=(:black,1,:+),aspect_ratio=:equal,xlims=zerone,ylims=zerone),
        		 scatter3d(randPnts[:,1],randPnts[:,2],res[:],title="max error: $(round(maximum(abs,res .- truth),digits = 5))",m=(:black,1,0.2),zlims=(0,10),camera=(30,70),xlims=zerone,ylims=zerone),
        		 scatter3d(spPoints[:,1],spPoints[:,2],pred[:],title="grid prediction",m=(:red,1,0.2),zlims=(0,10),camera=(30,70),xlims=zerone,ylims=zerone,zforeground_color_grid=:black),
        		 layout=(1,3),leg=false
        		 )
        end
        if save
            gif(anim,joinpath(dirname(@__FILE__),"ex3.gif"),fps=1)
        end
        return 0

    end