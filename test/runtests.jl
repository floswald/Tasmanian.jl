using Tasmanian
using Test

# write your own tests here
@testset "Basics" begin

	dim = 2
	out = 1
	depth = 5
	tsg = Tasmanian.TasmanianSG(dim,out,depth)
	@test tsg.version == VersionNumber("8.0")
	Tasmanian.makeLocalPolynomialGrid!(tsg)
	@test Tasmanian.isLocalPolynomial(tsg)
	@test !Tasmanian.isGlobal(tsg)
	@test Tasmanian.getDims(tsg) == Tasmanian.getNumDimensions(tsg)
	@test Tasmanian.getNout(tsg) == Tasmanian.getNumOutputs(tsg)
	@test size(Tasmanian.getPoints(tsg)) == (dim, 145)
end

@testset "run basic" begin
    t = Tasmanian.run()
    @test isa(t,TasmanianSG)
end

@testset "run examples" begin
    Tasmanian.ex1()
    Tasmanian.ex2()
    Tasmanian.ex3()
end

@testset "BasciIO" begin
    include("testBasicIO.jl")
    include("testCommon.jl")
    checkCopySubgrid()
end
