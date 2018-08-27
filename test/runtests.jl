using Tasmanian
using TestSetExtensions
using Test

# write your own tests here
@testset "Basics" begin

	dim = 2
	out = 1
	depth = 5
	tsg = Tasmanian.TasmanianSG(dim,out,depth)
	@test tsg.version == VersionNumber("5.1.0")
	Tasmanian.makeLocalPolynomialGrid!(tsg)
	@test Tasmanian.isLocalPolynomial(tsg)
	@test !Tasmanian.isGlobal(tsg)
	@test Tasmanian.getDims(tsg) == Tasmanian.getNumDimensions(tsg)
	@test Tasmanian.getNout(tsg) == Tasmanian.getNumOutputs(tsg)
	@test size(Tasmanian.getPoints(tsg)) == (145,dim)
end
