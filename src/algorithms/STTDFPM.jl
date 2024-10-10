# function STTDFPM(F::Function,P::Function, x0::Vector{Float64}, maxIter::Int64, μ1::Float64, μ2::Float64, λ::Float64, γ::Float64, η::Float64, ϕ::Float64, tol::Float64, ρ::Float64, ζ::Float64)
#     x = x0;
#     fx = F(x0);
#     normfx=norm(fx)
#     f_evals =1
#     d = -F(x0);
#     normd=norm(d);
#    iter = 0;
#     #z = x0  
#     while normfx > tol && iter< maxIter
#         α =0.000;
#         n=1000;

#         r2 = min(max(μ1, normfx), min(normfx, μ2))
#         for i in 1:n
#             fxd= F(x + ρ * η^i*d)
#             f_evals +=1
#             if -dot(fxd, d) ≥ ζ * ρ * η^i *r2* (normd)^2
#                 α = max( ρ * η^i);
#                 #println("α = ",  ρ * η^i)
#                 break
#             end

#         end

#         if  normfx <= tol
#             #println("Converged at x = ", x, " with F(x) = ", fx)
#             return normfx, iter, f_evals,x 
#             else
#                 z = x + α * d
#             #println("z = ", z)
#             fz=F(z);
#             normfz = norm(fz)
#             f_evals +=1
#                 if normfz <= tol
#                    # println("Converged at z = ", z, " with F(z) = ", F(z))
#                     return normfz, iter, f_evals, z
#                     else
#                         t = x -  γ*(dot(fz, x - z) / (normfz)^2) * fz;

#                 end

#         xnew = P(t);       
#         fxnew=F(xnew);
#         normfxnew=norm(fxnew);
#         f_evals +=1;
#         y = fxnew - fx; 
#         s = xnew - x + λ * y;
#         normy=norm(y);
#         θ = dot(s, y) / (normy)^2;
#         r1 = min(max(μ1, θ), min(θ, μ2));
#         v = max(ϕ*(normd)*(normy), (normfxnew)^2);
#         β = dot(fxnew,y) /v;
#         τ = dot(fxnew,d)/v;
#         d = -r1 *fxnew + β * d - τ*y;

#         x = xnew;
#         fx = F(x);
#         normfx=norm(fx)
#         f_evals +=1;
#         iter +=1;


#         end




#        # println(" x = x, F(x) = $fx")
#     end

#     #println("Maximum iterations reached. Final x = ", x, ", F(x) = ", F(x))
#     return normfx, iter, f_evals, x


# end

# function STTDFPM1(F::Function, x0::Vector{Float64}, maxIter::Int64, μ1::Float64, μ2::Float64, α1::Float64, α2::Float64, λ::Float64, γ::Float64, η::Float64, ϕ::Float64, tol::Float64, ρ::Float64, ζ::Float64)
#     x = x0;
#     fx = F(x0);
#     normfx=norm(fx)
#     f_evals =1
#     d = -F(x0);
#     normd=norm(d);
#    iter = 0;
#     #z = x0  
#     while normfx > tol && iter< maxIter
#         α =0.000;
#         n=1000;

#         r2 = min(max(μ1, normfx), min(normfx, μ2))
#         for i in 1:n
#             fxd= F(x + ρ * η^i*d)
#             f_evals +=1
#             if -dot(fxd, d) ≥ ζ * ρ * η^i *r2* (normd)^2
#                 α = max( ρ * η^i);
#                 #println("α = ",  ρ * η^i)
#                 break
#             end

#         end

#         if  normfx <= tol
#             #println("Converged at x = ", x, " with F(x) = ", fx)
#             return normfx, iter
#             break
#             else
#                 z = x + α * d
#             #println("z = ", z)
#             fz=F(z);
#             normfz = norm(fz)
#             f_evals +=1
#                 if normfz <= tol
#                    # println("Converged at z = ", z, " with F(z) = ", F(z))
#                     return normfz, iter, f_evals
#                     else
#                         t = x -  γ*(dot(fz, x - z) / (normfz)^2) * fz;

#                 end



#         end

#         xnew = P*t;       
#         fxnew=F(xnew);
#         normfxnew=norm(fxnew);
#         f_evals +=1;
#         y = fxnew - fx; 
#         s = xnew - x + λ * y;
#         normy=norm(y);
#         θ = dot(s, y) / (normy)^2;
#         r1 = min(max(α1, θ), min(θ, α2));
#         v = max(ϕ*(normd)*(normy), (normfxnew)^2);
#         β = dot(fxnew,y) /v;
#         τ = dot(fxnew,d)/v;
#         d = -r1 *fxnew + β * d - τ*y;

#         x = xnew;
#         fx = F(x);
#         normfx=norm(fx)
#         f_evals +=1;
#         iter +=1;



#        # println(" x = x, F(x) = $fx")
#     end

#     #println("Maximum iterations reached. Final x = ", x, ", F(x) = ", F(x))
#     return normfx, iter, f_evals


# end
"""
SOURCE: Ibrahim, A.H., Alshahrani, M. & Al-Homidan, S.
Two classes of spectral three-term derivative-free method for solving nonlinear 
equations with application. 
Numer Algor 96,1625–1645 (2024). 
https://doi.org/10.1007/s11075-023-01679-7

"""

# function STTDFPM(F::Function,
#     x0::Vector{Float64},
#     P::Function, Pcheck::Function,
#     maxIter::Int64=2_000,
#     tol::Float64=1e-11,
#     η1::Float64=0.001,
#     η2::Float64=0.6,
#     αmin::Float64=1e-10,
#     αmax::Float64=1e30,
#     r::Float64=0.1,
#     γ::Float64=1.8,
#     β::Float64=0.5,
#     ϕ::Float64=0.2,
#     σ::Float64=0.01
# )
#     x = x0
#     z = zeros(length(x))
#     fx = F(x0)
#     normfx = norm(fx)
#     f_evals = 1
#     d = -fx
#     normd = norm(d)
#     fz = F(z)
#     normfz = norm(fz)
#     f_evals += 1
#     iter = 0
#     τk = 0.0
#     n = 1_000


#     while true
#         if iter > maxIter
#             return normfx, iter, f_evals, x, :max_iters_reached
#         end
#         if normfx <= tol
#             return normfx, iter, f_evals, x, :converged_with_x
#         end

#         if normd <= tol
#             return normfx, iter, f_evals, x, :converged_when_d
#         end
#         r2 = min(max(η1, normfz), min(normfz, η2))
#         for i in 1:n
#             τk = β^i
#             fxd = F(x + τk * d)
#             f_evals += 1
#             if -dot(fxd, d) ≥ σ * τk * r2 * normd^2

#                 break
#             end
#         end

#         # Check for convergence
#         z = x + τk * d
#         fz = F(z)
#         normfz = norm(fz)
#         f_evals += 1

#         if Pcheck(z) && normfz <= tol
#             return normfz, iter, f_evals, z, :converged_with_z
#         end
#         μ = (dot(fz, -τk * d) / normfz^2)
#         t = x - γ * μ * fz


#         # Projection step
#         xnew = P(t)
#         fxnew = F(xnew)
#         #normfxnew = norm(fxnew)
#         f_evals += 1

#         y = fxnew - fx
#         s = xnew - x + r * y
#         normy = norm(y)
#         θ = dot(s, y) / normy^2
#         r1 = min(max(αmin, θ), min(θ, αmax))
#         v = max(ϕ * normd * normy, normfx^2)
#         βk = dot(fxnew, y) / v
#         τ = dot(fxnew, d) / v
#         d = -r1 * fxnew + βk * d - τ * y

#         # Update for next iteration
#         x = xnew
#         fx = F(x)
#         normfx = norm(fx)
#         f_evals += 1
#         iter += 1
#     end

#     # Return result after max iterations

# end

function STTDFPM(F::Function,
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
