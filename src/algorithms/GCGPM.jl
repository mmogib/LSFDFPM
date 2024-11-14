function GCGPM(F::Function, x0::Vector{Float64},
    P::Function,
    Pchcek::Function, args...)
    GCGPM(F, x0, P, args...)
end
function GCGPM(F::Function,
    x0::Vector{Float64},
    P::Function,
    maxIter::Int64=2_000,
    tol::Float64=1e-11,
    𝛷::Float64=1.48,
    σ2::Float64=0.3,
    β::Float64=1.0,
    ρ::Float64=0.6,
    σ::Float64=0.001
)
    x = x0
    fx = F(x0)
    normfx = norm(fx)
    f_evals = 1
    iter = 0
    d = -fx
    normd = norm(d)
    α = 0.000
    n1 = 1000
    while true
        if iter > maxIter
            return normfx, iter, f_evals, x, :max_iters_reached
        end
        if normfx <= tol
            return normfx, iter, f_evals, x, :converged_with_x
        end


        for i in 1:n1
            α = β * ρ^i
            fxd = F(x + α * d)
            normfxd = norm(fxd)
            f_evals += 1
            if -dot(fxd, d) ≥ σ * α * normfxd * normd^2
                break
            end
        end
        z = x + α * d
        fz = F(z)
        normfz = norm(fz)
        f_evals += 1
        τk = dot(fz, -α * d) / (normfz)^2
        t = x - τk * fz


        # Apply the projection matrix P to compute xnew
        xnew = P(t)
        fxnew = F(xnew)
        f_evals += 1
        y = fxnew - fx
        λ = 1 + max(0, -dot(y, d) / normd^2)
        w = y + λ * d
        k = dot(w, d)
        c = dot(fxnew, d) / k
        βk = (dot(fxnew, w) - 𝛷 * c * norm(w)^2) / k
        d = -𝛷 * fxnew + βk * d + σ2 * c * w
        x = xnew
        fx = F(x)
        normfx = norm(fx)
        normd = norm(d)
        iter += 1
    end

end