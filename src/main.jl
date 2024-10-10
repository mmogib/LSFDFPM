include("depts.jl")
include("includes.jl")
disp = vscodedisplay

x0 = zeros(100) .+ 0.1
x0[1] = -2
ψ(x) = x .^ 2 .- 1
μ0 = 0.5
α = 0.54
c1 = c2 = 0
γ = 0.001
maxitrs = 2_000
ε = 1e-12
x = LSFDFPM(x0, ψ, μ0, α, ε, maxitrs, c1, c2, γ)

norm(F4(x))
x