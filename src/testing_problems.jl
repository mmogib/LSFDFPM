include("depts.jl")
include("initial_points.jl")
include("problems.jl")
using BenchmarkTools

f1 = F1
f2 = getfield(Main, Symbol("$(String(Symbol(f1)))v"))
dm = 1000
points = map(x -> x[3], values(getInitialPoints(dm)))

@benchmark f1.(points)
@benchmark f2.(points)
if norm(map(norm, f1.(points) - f2.(points)), Inf) ≈ 0
    println("same")
else
    print("different")
end

testThem(F1, F1v, 1000)