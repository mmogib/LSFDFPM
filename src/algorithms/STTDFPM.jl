"""
SOURCE: Ibrahim, A.H., Alshahrani, M. & Al-Homidan, S.
Two classes of spectral three-term derivative-free method for solving nonlinear 
equations with application. 
Numer Algor 96,1625–1645 (2024). 
https://doi.org/10.1007/s11075-023-01679-7

"""
function STTDFPM(
    F::Function,
    x0::Vector{Float64},
    P::Function, Pcheck::Function,
    maxiters::Int64=2_000,
    ϵ::Float64=1e-11,
    η::Float64=0.001,
    ξ::Float64=0.6,
    αmin::Float64=1e-10,
    αmax::Float64=1e30,
    r::Float64=0.1,
    γ::Float64=1.8,
    β::Float64=0.5,
    ψ::Float64=0.2,
    σ::Float64=0.01)
    # t, β, σ, γ, η, ξ = params.t, params.β, params.σ, params.γ, params.η, params.ξ
    # αmin, αmax, r, ψ = params.αmin, params.αmax, params.r, params.ψ
    # ϵ, maxiters, stop_criterion = options.tol, options.maxiters, options.stopping
    evals = 0
    x1 = x0
    Fx1 = F(x1)
    Fx1_norm = norm(Fx1)
    evals += 1
    d1 = -Fx1
    k = 0

    while true
        if Fx1_norm <= ϵ
            return Fx1_norm, k, evals, x1, :converged_with_x
        end
        if k > maxiters
            return Fx1_norm, k, evals, x1, :max_iters_reached
        end
        τ1, tevals = findtk(F, d1, x1, σ, β, η, ξ)
        evals = evals + tevals
        if isnothing(τ1)
            return Fx1_norm, k, evals, x1, :failed_line_search
        end

        z1 = x1 + τ1 * d1
        Fz1 = F(z1)
        evals = evals + 1
        Fz1_norm = norm(Fz1)
        if Pcheck(z1) && Fz1_norm <= ϵ
            return Fz1_norm, k, evals, z1, :converged_with_z
        end
        μ1 = dot(Fz1, x1 - z1) / (Fz1_norm^2)
        x1, x0, Fx0 = P(x1 - γ * μ1 * Fz1), x1, Fx1

        Fx1 = F(x1)
        evals += 1
        Fx1_norm = norm(Fx1)
        d1 = findDk(Fx0, Fx1, d1, x0, x1, r, ψ, αmin, αmax)
        if 10 * norm(d1) <= ϵ
            return Fx1_norm, k, evals, x1, :converged_with_x
        end
        k = k + 1
    end


end

function findDk(Fu0, Fu1, d0, u0, u1, r, ψ, αmin, αmax)
    y0 = Fu1 - Fu0
    s0 = u1 - u0 + r * y0
    # v1 = max(ψ * norm(d0) * norm(y0), dot(d0, y0), norm(Fu0)^2)
    v1 = max(ψ * norm(d0) * norm(y0), norm(Fu0)^2)
    β1 = dot(Fu1, y0) / v1
    α12 = dot(Fu1, d0) / v1
    α11 = min(αmax, max(αmin, dot(s0, y0) / dot(y0, y0)))
    return -α11 * Fu1 + β1 * d0 - α12 * y0
end


function findtk(F::Function, d::Vector{Float64}, u::Vector{Float64}, σ::Float64, β::Float64, η::Float64, ξ::Float64)
    # ϵ = 1e-6
    max_iters = 10_000
    for i in 0:max_iters
        tk = β^i
        arg = F(u + tk .* d)
        lhs = -dot(arg, d)
        χ = norm(arg)
        Pηξχ = min(ξ, max(χ, η))
        rhs = σ * tk * Pηξχ * norm(d)^2
        # rhs = σ * tk * norm(d)^2
        if lhs >= rhs || tk <= 1e-5
            return tk, i
        end
    end
    nothing, max_iters
end
