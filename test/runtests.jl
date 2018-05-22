using Tasmanian
using TestSetExtensions
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
@testset "Basics" begin
	tsg = Tasmanian.TasmanianSG()
	@test tsg.version == VersionNumber("6.0.0")
	dim = 2
	out = 1
	depth = 1
	Tasmanian.makeLocalPolynomialGrid!(tsg,dim,out,depth)
	@test Tasmanian.getNumPoints(tsg) == out
	@test Tasmanian.getNumDimensions(tsg) == dim
	@test Tasmanian.getNumLoaded(tsg) == 0
	@test Tasmanian.getPoints(tsg) == [0.0 0.0]
end
