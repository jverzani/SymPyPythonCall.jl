using SymPyPythonCall
import Symbolics
using Test


@testset "Symbolics" begin
    @syms x
    𝑥 = SymPyPythonCall.PythonCall.pyconvert(Symbolics.Num, x)
    @test isa(𝑥, Symbolics.Num)
end
