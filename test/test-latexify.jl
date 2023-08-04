using SymPyPythonCall
import SymPyPythonCall: Latexify
using Test

@testset "Latexify" begin
  @syms α
  @test_nowarn Latexify.latexify(- 3/(8*α) + 1/(8*α^2))
  @test  occursin("\\cdot", Latexify.latexify(- 3/(8*α) + 1/(8*α^2), cdot=true))
  @test !occursin("\\cdot", Latexify.latexify(- 3/(8*α) + 1/(8*α^2), cdot=false))
end
