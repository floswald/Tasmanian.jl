using Tasmanian
using TestSetExtensions
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
@testset "Basics" begin

	dim = 2
	out = 1
	depth = 5
	tsg = Tasmanian.TasmanianSG(dim,out,depth)
	@test tsg.version == VersionNumber("6.0.0")
	Tasmanian.makeLocalPolynomialGrid!(tsg)
	@test Tasmanian.isLocalPolynomial(tsg)
	@test !Tasmanian.isGlobal(tsg)
	@test Tasmanian.getDims(tsg) == Tasmanian.getNumDimensions(tsg)
	@test Tasmanian.getNout(tsg) == Tasmanian.getNumOutputs(tsg)
	@test size(Tasmanian.getPoints(tsg)) == (145,dim)
end
