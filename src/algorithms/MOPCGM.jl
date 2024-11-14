"""
SOURCE: Sabi’u, J., Shah, A., Stanimirović, P. S., Ivanov, B., & Waziri, M. Y. 
(2023). 
Modified optimal Perry conjugate gradient method for solving system of 
    monotone equations with applications. Applied Numerical Mathematics, 
    184, 431–445. 
    https://doi.org/10.1016/j.apnum.2022.10.016
"""
function MOPCGM(F::Function, x0::Vector{Float64},
    P::Function, Pcheck::Function,
    maxIter::Int64=2_000,
    tol::Float64=1e-11,
    λ::Float64=0.1,
    η::Float64=0.9,
    ρ::Float64=0.1,
    ζ::Float64=0.0001
)
    x = x0
    n1 = 1_000
    fx = F(x0)
    f_evals = 1
    d = -F(x0)
    iter = 0
    normfx = norm(fx)
    normd = norm(d)
    α = 0
    while true
        if iter > maxIter
            return normfx, iter, f_evals, x, :max_iters_reached
        end
        if normfx <= tol
            return normfx, iter, f_evals, x, :converged_with_x
        end

        for i in 1:n1
            α = ρ * η^i
            fxd = F(x + α * d)
            f_evals += 1
            if -dot(fxd, d) ≥ ζ * α * (normd)^2
                break
            end
        end

        z = x + α * d
        fz = F(z)
        normfz = norm(fz)
        f_evals += 1

        if Pcheck(z) && normfz <= tol
            return normfz, iter, f_evals, z, :converged_with_z
        end
        bk = dot(fz, -α * d) / (normfz)^2
        t = x - bk * fz
        xnew = P(t)
        s = z - x
        fxnew = F(xnew)
        f_evals += 1
        normfxnew = norm(fxnew)
        y = fxnew - fx + λ * s

        θ = dot(s, y) / norm(s)^2
        β = dot(y - θ * s, fxnew) / dot(d, y)

        r = dot(fxnew, d) / (normfxnew)^2
        d = -(1 + β * r) * fxnew + β * d
        x = xnew
        fx = F(x)
        f_evals += 1
        normd = norm(d)
        normfx = norm(fx)
        iter += 1
    end
end

#function MOPCGM2(F::Function,P::Function, maxIter::Int64, x0::Vector{Float64}, 𝛷::Float64, λ::Float64, η::Float64, tol::Float64, ρ::Float64, ζ::Float64)
#   x = x0;
#   fx = F(x0);
#    f_evals=1;
#    d = -F(x0);
#    iter = 0;
#    normfx = norm(fx);

#    while normfx > tol && iter < maxIter
#        α =0.000;
#        n1=1000;
#        for i in 1:n1
#            fxd=F(x + ρ * η^i*d);
#            f_evals +=1;
#            if -dot(fxd, d) ≥ ζ * ρ * η^i * (normfx)^2
#                α = max( ρ * η^i);
#                #println("α = ",  ρ * η^i)
#                break
#            end

#        end

#        if  normfx <= tol
#println("Converged at x = ", x, " with F(x) = ", fx)
#            return normfx, iter, f_evals
#            else
#                z = x + α * d;
#println("z = ", z)
#            fz=F(z);
#            normfz=norm(fz);
#            f_evals +=1;

#                if normfz <= tol 
#println("Converged at z = ", z, " with F(z) = ", F(z))
#                    return normfz, iter, f_evals
#                    else
#                        t = x - (dot(fz, x - z) / (normfz)^2) * fz;

#                end

#        xnew = P(t);      
#       s = xnew - x;
#       fxnew =F(xnew);
#        normfxnew = norm(fxnew);
#        f_evals +=1;
#        y = fxnew - fx + λ * s;
#        r = dot(fx, d)/(normfxnew)^2;
#        θ = 𝛷* dot(s, y) / norm(s)^2;
#        β = dot(y - θ * s, fxnew) / dot(d, y);
#        d =  -𝛷* fxnew - β * r * fxnew + β * d;
#        x = xnew;
#        fx = F(x);
#        normfx=norm(fx);
#        f_evals +=1
#        iter +=1



#        end




# println(" x = , F(x) = $fx")
#   end

#println("Maximum iterations reached. Final x = ", x, ", F(x) = ", F(x))
#    return norm(F(x)), iter, f_evals


#end

