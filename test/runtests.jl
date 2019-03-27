using Tasmanian
using Test

# write your own tests here
@testset "Basics" begin

	dim = 2
	out = 1
	depth = 5
	tsg = Tasmanian.TasmanianSG(dim,out,depth)
	@test tsg.version == VersionNumber("6.0")
	Tasmanian.makeLocalPolynomialGrid!(tsg)
	@test Tasmanian.isLocalPolynomial(tsg)
	@test !Tasmanian.isGlobal(tsg)
	@test Tasmanian.getDims(tsg) == Tasmanian.getNumDimensions(tsg)
	@test Tasmanian.getNout(tsg) == Tasmanian.getNumOutputs(tsg)
	@test size(Tasmanian.getPoints(tsg)) == (145,dim)
end

@testset "run basic" begin
    t = Tasmanian.run()
    @test isa(t,TasmanianSG)
end

@testset "run examples" begin
    @test Tasmanian.ex1() == 0
    @test Tasmanian.ex2() == 0
    @test Tasmanian.ex3() == 0
end
