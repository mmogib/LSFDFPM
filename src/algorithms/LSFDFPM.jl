function P_H(v::Vector{<:Real}, x::Vector{<:Real})
    v_norm2 = dot(v, v)
    return (u::Vector{<:Real}) -> u + (dot(-v, u - x) / v_norm2)v
end


"""
LSFDFPM(
    ψ::Function,
    x0::Vector{Float64},
    SD::Function,
    LS::Union{Nothing,Function},
    maxiters::Int64=2_000,
    ϵ::Float64=1e-11,
    time_limit::Union{Nothing,Integer}=nothing,
    μ0::Float64=0.5,
    γ::Float64=1.8,
    α::Union{Nothing,Float64}=nothing,
)


SOURCE: Ibrahim, A.H., Alshahrani, M. & Al-Homidan, S.
Two classes of spectral three-term derivative-free method for solving nonlinear 
equations with application. 
Numer Algor 96,1625–1645 (2024). 
https://doi.org/10.1007/s11075-023-01679-7

"""
function LSFDFPM(
    ψ::Function,
    x0::Vector{Float64},
    SD::Function,
    LS::Union{Nothing,Function},
    maxiters::Int64=2_000,
    ϵ::Float64=1e-11,
    time_limit::Union{Nothing,Integer}=nothing,
    μ0::Float64=0.5,
    γ::Float64=1.8,
    α::Union{Nothing,Float64}=nothing,
)
    # gamma = 0.2
    # alpha = 0.5
    # r = 0.1
    # psi = 0.2
    # mu0 = 0.1
    # mu = mu0
    Fe = 0
    x1 = copy(x0)

    ψx1 = ψ(x1)
    ψx1_norm = norm(ψx1)
    Fe += 1
    d1 = -ψx1
    k = 0
    μ1 = μ0
    # line search
    lsparams = LSParams(1, Inf64, Inf64, Inf64)
    Random.seed!(2024)
    N_d = truncated(Normal(), 0.2, 0.9999)
    ls = if isnothing(LS)
        α0 = isnothing(α) ? rand(N_d) : α
        (_x, _y; params::LSParams) -> LSOutput(α0, 0, Inf64, Inf64)
    else
        LS
    end
    time_start = Dates.now()
    time_criterion() = begin
        if isnothing(time_limit)
            return false
        else
            current_time = Dates.now()
            if round(current_time - time_start, Second) > Second(time_limit)
                return true
            end
            return false
        end
    end

    while true
        if ψx1_norm <= ϵ
            return ψx1_norm, k, Fe, :converged_with_x
        end
        if k > maxiters
            return ψx1_norm, maxiters, Fe, :max_iters_reached
        end
        if time_criterion()
            return ψx1_norm, maxiters, Fe, :max_time_reached
        end
        lsout = ls(d1, x1; params=lsparams)
        τ1, tFe, Δ, αstab = lsout.α, lsout.iters, lsout.Δ, lsout.αstab

        Fe = Fe + tFe
        if isnothing(τ1)
            return ψx1_norm, k, Fe, :failed_line_search
        end

        z1 = x1 + τ1 * d1
        ψz1 = ψ(z1)
        Fe = Fe + 1
        ψz1_norm = norm(ψz1)
        if ψz1_norm <= ϵ
            return ψz1_norm, k, Fe, :converged_with_z
        end

        if isnothing(LS)
            u = x1 - μ1 * ψz1  # A placeholder for the projection

            v = x1 - μ1 * ψx1 - z1
            x_next = P_H(v, z1)(u)

            # Step 5: Update μ_j+1 based on condition
            norm_ψ_x_ψ_z = norm(ψx1 - ψz1)
            if norm_ψ_x_ψ_z > 0
                μ1 = min(γ * norm(τ1 * d1) / norm_ψ_x_ψ_z, μ1)
            end

            # Update variables for the next iteration
            x1, x0, ψx0 = x_next, x1, ψx1
        else

            μ1 = dot(ψz1, x1 - z1) / (ψz1_norm^2)
            x1, x0, ψx0 = x1 - γ * μ1 * ψz1, x1, ψx1
        end
        ψx1 = ψ(x1)
        Fe += 1
        ψx1_norm = norm(ψx1)
        d1 = SD(ψx0, ψx1, d1, x0, x1, k)
        # if 10 * norm(d1) <= ϵ
        #     return ψx1_norm, k, Fe, :converged_with_x
        # end
        k = k + 1
        lsparams = LSParams(k, Δ, 0.0, αstab)
    end


end

# function findDk(Fu0, Fu1, d0, u0, u1, r, ψ, αmin, αmax)
#     y0 = Fu1 - Fu0
#     s0 = u1 - u0 + r * y0
#     # v1 = max(ψ * norm(d0) * norm(y0), dot(d0, y0), norm(Fu0)^2)
#     v1 = max(ψ * norm(d0) * norm(y0), norm(Fu0)^2)
#     β1 = dot(Fu1, y0) / v1
#     α12 = dot(Fu1, d0) / v1
#     α11 = min(αmax, max(αmin, dot(s0, y0) / dot(y0, y0)))
#     return -α11 * Fu1 + β1 * d0 - α12 * y0
# end


# function findtk(ψ::Function, d::Vector{Float64}, u::Vector{Float64}, σ::Float64, β::Float64, η::Float64, ξ::Float64)
#     # ϵ = 1e-6
#     max_iters = 10_000
#     for i in 0:max_iters
#         tk = β^i
#         arg = ψ(u + tk .* d)
#         lhs = -dot(arg, d)
#         χ = norm(arg)
#         Pηξχ = min(ξ, max(χ, η))
#         rhs = σ * tk * Pηξχ * norm(d)^2
#         # rhs = σ * tk * norm(d)^2
#         if lhs >= rhs || tk <= 1e-5
#             return tk, i
#         end
#     end
#     nothing, max_iters
# end
