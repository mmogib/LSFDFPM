# """
# Source:
# Problem 4.6 in
# Jamilu Sabi'u, Abdullah Shah, Predrag S. Stanimirović, Branislav Ivanov, Mohammed Yusuf Waziri, 
# Modified optimal Perry conjugate gradient method for solving system of monotone equations with applications,Applied Numerical Mathematics,
# Volume 184, 2023, Pages 431-445, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2022.10.016.

# Projected on [0, ∞] or [-1, ∞]

# name: PolynomialSineCosine
# """
# function F1(x)
#     N = 1 / length(x)
#     s = N * sum(x)
#     x .* cos.(x .- N) .* (
#         sin.(x) .- 1 .- (x .- 1) .^ 2 .- s
#     )
# end

# """
# Source: 
# Problem 5.2
# Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
# Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96. https://www.global-sci.org/intro/article_detail/jcm/8675.html#
# Example: 4.1 in
# Wang, C., Wang, Y. & Xu, C. A projection method for a system of nonlinear monotone equations with convex constraints. 
# Math Meth Oper Res 66, 33–46 (2007). https://doi.org/10.1007/s00186-006-0140-y

# name: ExponetialI
# """
# function F2(x)
#     exp(x) - 1
# end



# """
# Source:
# Problem 3 in Appendix
# William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
# Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493

# name: ExponetialIII
# """
# function F3(x)
#     X2 = x .^ 2
#     Ii = [i / 10 for i in 1:length(x)]
#     Y = vcat(
#         X2[1:end-1] .+ exp.(-X2[1:end-1]),
#         exp(-X2[end])
#     )
#     Ii .* (1 .- Y)
# end


# """
# Source:
# Problem 9 
# J. Sabi’u, A. Shah, M.Y. Waziri, M.K. Dauda, A new hybrid approach for solving large-scale monotone nonlinear equations, 
# J. Math. Fund. Sci. 52 (2020) 17–26 https://doi.org/10.5614/j.math.fund.sci.2020.52.1.2

# name: PolynomialI
# """
# function F4(x)
#     vcat(
#         4 .* x[1:end-1] .+ (x[2:end] .- 2 .* x[1:end-1]) .- ((x[2:end] .^ 2) / 3),
#         4 * x[end] + (x[end-1] - 2 * x[end]) - x[end-1]^2 / 3
#     )
# end



# """
# Source:
# Problem 2 in
# Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
# Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.

# Project to Rn₊

# Modified Problem 1. in
# W.J. Zhou, D.H. Li, A globally convergent BFGS method for nonlinear monotone equations without any merit functions, 
# Math. Comput. 77 (264) (2008) 2231–2240.
# Projected on [-2, ∞]

# name: SmoothSine
# """
# function F5(x)
#     v = try

#         2 * x - sin(x)
#     catch
#         println(max(x...), min(x...))
#         throw("problem with this z")
#     end
#     v
# end



# """
# Source:
# Problem 1 in
# Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
# Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.tional and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

# Project to Rn₊

# name: NonsmoothSine
# """
# function F6(x)
#     2 * x - sin(abs(x))
# end


# """
# Source:
# Problem 2 in
# Zhensheng Yu, Ji Lin, Jing Sun, Yunhai Xiao, Liying Liu, Zhanhui Li, Spectral gradient projection method for monotone nonlinear equations with convex constraints,
# Applied Numerical Mathematics, Volume 59, Issue 10, 2009, Pages 2416-2423, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2009.04.004.

# with Projection on 
#     C = {x ∈ Rⁿ | ∑xᵢ ≤ n, xᵢ≥ -1}

# Modified from Problem 1 in
# Li Zhang, Weijun Zhou, Spectral gradient projection method for solving nonlinear monotone equations, 
# Journal of Computational and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

# name: ModifiedNonsmoothSine
# """
# function F7(x)
#     x - sin(abs(x - 1))
# end

# """
# Source:
# Problem 2
# Gao PT, He CJ (2018) An efficient three-term conjugate gradient method for nonlinear monotone equations
# with convex constraints. Calcolo. https://doi.org/10.1007/s10092-018-0291-2


# with Projection on 
#     C = {x ∈ Rⁿ | ∑xᵢ ≤ n, xᵢ≥ -1}

# Modified from Problem 1 in
# Li Zhang, Weijun Zhou, Spectral gradient projection method for solving nonlinear monotone equations, 
# Journal of Computational and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

# name: ModifiedNonsmoothSine2
# """
# function F8(x)
#     x - sin(abs(x) - 1)
# end


# """
# Source:
# Problem 1
# Gao, P., He, C. An efficient three-term conjugate gradient method for nonlinear monotone equations with convex constraints. 
# Calcolo 55, 53 (2018). https://doi.org/10.1007/s10092-018-0291-2

# with Projection on 
#     Rn+

# name: ExponetialSineCosine
# """
# function F9(x)
#     (exp.(x)) .^ 2 + 3 * sin.(x) * cos.(x) - 1
# end

# """
# Source: 
# Problem 4.2 in
# Li Zheng, Lei Yang, Yong Liang, A conjugate gradient projection method for solving equations with convex constraints,
# Journal of Computational and Applied Mathematics, Volume 375, 2020, 112781, ISSN 0377-0427,
# https://doi.org/10.1016/j.cam.2020.112781.
# Modified from Problem 3
# Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
# Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.

# name: ModifiedTrigI
# """
# function F10(x)
#     vcat(
#         x[1] + sin(x[1]) - 1,
#         -x[1:end-2] + 2 * x[2:end-1] + sin.(x[2:end-1]) .- 1,
#         x[end] + sin(x[end]) - 1
#     )
# end

# """
# Source:
# Problem 10 in 
# Y. Bing, G. Lin, An efficient implementation of Merrill’s method for sparse or partially separable systems of nonlinear equations, 
# SIAM. J. Optim.1 (2) (1991) 206–221.

# name: Tridiagonal
# """
# function F11(x)
#     # Ii = [i for i in 1:length(x)]
#     x0 = vcat(0, x)
#     x1 = vcat(x, 0)
#     h = 1 / (1 + length(x))
#     x .- exp.(cos.(h .* (x0[1:end-1] .+ x .+ x1[2:end])))
# end


# """
# Source:
# Problem 4.4 in
# Li Zheng, Lei Yang, Yong Liang, A conjugate gradient projection method for solving equations with convex constraints,
# Journal of Computational and Applied Mathematics, Volume 375, 2020, 112781, ISSN 0377-0427,
# https://doi.org/10.1016/j.cam.2020.112781.
# Modified from Problem 10 in 
# Y. Bing, G. Lin, An efficient implementation of Merrill’s method for sparse or partially separable systems of nonlinear equations, 
# SIAM. J. Optim.1 (2) (1991) 206–221.

# name: ModifiedTridiagonal
# """
# function F12(x)
#     n = length(x)
#     vcat(
#         x[1] - exp(cos((x[1] + x[2]) / (n + 1))),
#         -x[2:end-1] - exp.(cos.((x[1:end-2] .+ x[2:end-1] .+ x[3:end]) ./ (n + 1))),
#         x[end] - exp(cos((x[end-1] + x[end]) / (n + 1)))
#     )
# end


# """
# Source:
# Problem 10 in Appendix
# William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
# Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493


# Problem 4.5 in
# Li Zheng, Lei Yang, Yong Liang, A conjugate gradient projection method for solving equations with convex constraints,
# Journal of Computational and Applied Mathematics, Volume 375, 2020, 112781, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2020.112781

# Projected on [-1, ∞]

# name: Logarithmic
# """
# function F13(x)
#     log.(x .+ 1) .- (x / length(x))
# end


# """
# Source:
# Problem 10 in Appendix
# William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
# Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493


# Modified in Problem 4.2 in
# Jamilu Sabi'u, Abdullah Shah, Predrag S. Stanimirović, Branislav Ivanov, Mohammed Yusuf Waziri, 
# Modified optimal Perry conjugate gradient method for solving system of monotone equations with applications,Applied Numerical Mathematics,
# Volume 184, 2023, Pages 431-445, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2022.10.016.

# Projected on [0, ∞] or [-1, ∞]

# name: NonmoothLogarithmic
# """
# function F14(x)
#     log.(abs.(x) .+ 1) .- (x / length(x))
# end



# # see https://www.cuter.rl.ac.uk//Problems/classification.shtml for classification

# # NAME: ARWHEAD         see (https://bitbucket.org/optrove/sif/raw/HEAD/ARWHEAD.SIF)
# # see also (https://vanderbei.princeton.edu/ampl/nlmodels/cute/arwhead.mod)
# #  classification OUR2-AN-V-0
# # O: none of the aboves
# # U: the problem is unconstrained
# # R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# # 2: the degree of the highest derivatives provided analytically within the problem description.
# # A: academic problem
# # N: the problem description does not contain any explicit internal variables
# # V: the number of variables in the problem can be chosen by the user
# # 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
# function ARWHEADFun(x)
#     # sum {i in 1..N-1} (-4*x[i]+3.0) + sum {i in 1..N-1} (x[i]^2+x[N]^2)^2
#     n = length(x)
#     sum((-4 * x[i] + 3.0) for i in 1:n-1) + sum((x[i]^2 + x[n]^2)^2 for i in 1:n-1)
# end
# # name: ARWHEADGrad
# function F15(x)
#     vcat(-4 .+ 4 * x[1:end-1] .* (x[1:end-1] .^ 2 .+ x[end]^2),
#         4 * x[end] * sum(x[1:end-1] .^ 2 .+ x[end]^2))
# end


# # name: PENALTY1 see (https://bitbucket.org/optrove/sif/raw/HEAD/PENALTY1.SIF) 
# # see also (https://vanderbei.princeton.edu/ampl/nlmodels/cute/penalty1.mod)
# # classification SUR2-AN-V-0
# # S: the objective function is a sum of squares
# # U: the problem is unconstrained
# # R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# # 2: the degree of the highest derivatives provided analytically within the problem description.
# # A: the problem is academic, that is, has been constructed specifically by researchers to test one or more algorithms,
# # N: the problem description does not contain any explicit internal variables.
# # V: the number of variables in the problem can be chosen by the user
# # 0: a nonnegative integer giving the actual (fixed) number of problem constraints.

# function PENALTY1Fun(x)
#     a = 1e-5
#     N = length(x)
#     # sum {i in 1..N} a*(x[i]-1)^2 + ( sum {j in 1..N} x[j]^2 - 1/4 )^2;
#     sum(a * (x[i] - 1)^2 for i in 1:N) + (sum(x[i]^2 for i in 1:N) - 0.25)^2
# end
# # name: PENALTY1Grad
# function F16(x)
#     t = sum(x .^ 2)
#     c = 1e-5
#     2 * c .* (x .- 1) + 4 * (t - 0.25) .* x
# end

# # NAME: DIXON3DQ see (https://bitbucket.org/optrove/sif/raw/HEAD/DIXON3DQ.SIF)
# # SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/dixon3dq.mod)
# # classification QUR2-AN-V-0
# # Q: the objective function is quadratic,
# # U: the problem is unconstrained
# # R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# # 2: the degree of the highest derivatives provided analytically within the problem description.
# # A: academic problem
# # N: the problem description does not contain any explicit internal variables
# # V: the number of variables in the problem can be chosen by the user
# # 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
# function DIXON3DQFun(x)
#     # (x[1]-1.0)^2 + sum {j in 2..n-1} (x[j]-x[j+1])^2 + (x[n]-1.0)^2
#     n = length(x)
#     (x[1] - 1.0)^2 + sum((x[j] - x[j+1])^2 for j in 2:n-1) + (x[n] - 1.0)^2
# end
# #name: DIXON3DQGrad
# function F17(x)
#     n = length(x)
#     return vcat(2 * (x[1] - 1),
#         2 * (x[2] - x[3]),
#         2 * (2 * x[3:n-1] - x[2:n-2] - x[4:n]),
#         2 * (2 * x[n] - x[n-1] - 1)
#     )
# end

# # NAME: GENHUMPS, see (https://bitbucket.org/optrove/sif/raw/HEAD/GENHUMPS.SIF)
# # SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/genhumps.mod)
# # classification OUR2-AN-V-0
# # O: none of the aboves
# # U: the problem is unconstrained
# # R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# # 2: the degree of the highest derivatives provided analytically within the problem description.
# # A: academic problem
# # N: the problem description does not contain any explicit internal variables
# # V: the number of variables in the problem can be chosen by the user
# # 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
# function GENHUMPSFun(x; ζ=20)
#     # sum {i in 1..N-1} ( sin (zeta*x[i])^2*sin(zeta*x[i+1])^2+0.05*(x[i]^2+x[i+1]^2) );
#     N = length(x)
#     sum((sin(ζ * x[i]) * sin(ζ * x[i+1]))^2 + 0.05 * (x[i]^2 + x[i+1]^2) for i in 1:(N-1))
# end
# # name: GENHUMPSGrad
# function F18(x; ζ=20)
#     TNTHX = 0.1 * x
#     ZX = ζ * x
#     ZX2 = 2 * ZX
#     SZX = (sin.(ZX)) .^ 2
#     SZX2 = sin.(ZX2)
#     ZSZX2 = ζ * SZX2
#     vcat(
#         ZSZX2[1] * SZX[2] + TNTHX[1],
#         ZSZX2[2:end-1] .* (SZX[1:end-2] .+ SZX[3:end]) .+ (2 * TNTHX[2:end-1]),
#         ZSZX2[end] * SZX[end-1] + TNTHX[end],
#     )
# end

# # NAME: ENGVAL1 see (https://bitbucket.org/optrove/sif/raw/HEAD/ENGVAL1.SIF)
# # SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/engval1.mod)
# # classification OUR2-AN-V-0
# # O: none of the aboves
# # U: the problem is unconstrained
# # R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# # 2: the degree of the highest derivatives provided analytically within the problem description.
# # A: academic problem
# # N: the problem description does not contain any explicit internal variables
# # V: the number of variables in the problem can be chosen by the user
# # 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
# function ENGVAL1Fun(x)
#     # sum {i in 1..N-1} (x[i]^2+x[i+1]^2)^2 + sum {i in 1..N-1} (-4*x[i]+3.0);
#     N = length(x)
#     sum((x[i]^2 + x[i+1]^2)^2 for i in 1:(N-1)) + sum(-4 * x[i] + 3.0 for i in 1:(N-1))
# end
# # name: ENGVAL1Grad
# function F19(x)
#     X2 = x .^ 2
#     vcat(
#         4 * (x[1] * (X2[1] + X2[2]) - 1),
#         4 * (x[2:end-1] .* (X2[1:end-2] + 2 * X2[2:end-1] + X2[3:end]) .- 1),
#         4 * x[end] * (X2[end-1] + X2[end])
#     )
# end
# # NAME: DIXMAANH see (https://bitbucket.org/optrove/sif/raw/HEAD/DIXMAANH.SIF)
# # SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/dixmaanh.mod)
# # classification OUR2-AN-V-0
# # O: none of the aboves
# # U: the problem is unconstrained
# # R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# # 2: the degree of the highest derivatives provided analytically within the problem description.
# # A: academic problem
# # N: the problem description does not contain any explicit internal variables
# # V: the number of variables in the problem can be chosen by the user
# # 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
# function DIXMAANHFun(x; α=1.0, β=0.26, γ=0.26, δ=0.26, K=[1 0 0 1])
#     N = length(x)
#     if N % 3 != 0
#         throw("The length of `x0` must be divisible by 3")
#     end
#     M = fld(N, 3)
#     1.0 + α * sum(x[i]^2 * (i / N)^K[1] for i in 1:N) +
#     β * sum(x[i]^2 * (x[i+1] + x[i+1]^2)^2 for i in 1:N-1) +
#     γ * sum(x[i]^2 * x[i+M]^4 for i in 1:2*M) +
#     δ * sum(x[i] * x[i+2*M] * (i / N)^K[4] for i in 1:M)
# end
# # name: DIXMAANHGrad
# function F20(x; α=1.0, β=0.26, γ=0.26, δ=0.26, K=[1 0 0 1])
#     N = length(x)
#     if N % 3 != 0
#         throw("The length of `x0` must be divisible by 3")
#     end
#     M = fld(N, 3)
#     X2 = x .^ 2
#     X4 = X2 .* X2
#     IOverN = [(i / N)^K[1] for i in 1:N]
#     IOverM = [(i / N)^K[4] for i in 1:M]
#     TwoAlphaXOverN = 2 * α * IOverN .* x
#     OnePlusX = 1 .+ x
#     OnePlus2X = OnePlusX .+ x
#     OnePlusX2 = OnePlusX .^ 2
#     G1 = TwoAlphaXOverN[1] + 2 * β * x[1] * X2[2] * OnePlusX2[2] + 2 * γ * x[1] * X4[1+M] + δ * IOverM[1] * x[1+2*M]
#     G2M = TwoAlphaXOverN[2:M] .+ 2 * β * (X2[1:M-1] .* x[2:M] .* OnePlusX[2:M] .* OnePlus2X[2:M] .+ x[2:M] .* X2[3:M+1] .* OnePlusX2[3:M+1]) .+
#           2 * γ * x[2:M] .* X4[2+M:2*M] .+ δ * IOverM[2:M] .* x[2+2*M:3*M]

#     # check
#     GM12M = TwoAlphaXOverN[M+1:2*M] .+
#             2 * β * (
#                 X2[M:2*M-1] .* x[M+1:2*M] .* OnePlusX[M+1:2*M] .* OnePlus2X[M+1:2*M] .+ x[M+1:2*M] .* X2[M+2:2*M+1] .* OnePlusX2[M+2:2*M+1]
#             ) .+
#             γ * (2 * x[M+1:2*M] .* X4[2*M+1:3*M] .+ 4 * X2[1:M] .* x[M+1:2*M] .* X2[M+1:2*M])


#     # check
#     G2M1Nm1 = TwoAlphaXOverN[2*M+1:N-1] .+
#               2 * β * (
#                   X2[2*M:N-2] .* x[2*M+1:N-1] .* OnePlusX[2*M+1:N-1] .* OnePlus2X[2*M+1:N-1] .+ x[2*M+1:N-1] .* X2[2*M+2:N] .* OnePlusX2[2*M+2:N]
#               ) .+
#               (γ * 4 * X2[M+1:2*M-1] .* x[2*M+1:N-1] .* X2[2*M+1:N-1]) .+
#               δ * IOverM[1:M-1] .* x[1:M-1]

#     GN = TwoAlphaXOverN[N] + 2 * β * (X2[N-1] * x[N] * OnePlusX[N] * OnePlus2X[N]) +
#          γ * (4 * X2[2*M] * x[N] * X2[N]) + δ * IOverM[M] .* x[M]
#     vcat(
#         G1,
#         G2M,
#         GM12M,
#         G2M1Nm1,
#         GN
#     )

# end

# function DIXMAANIFun(x)
#     DIXMAANHFun(x; α=1.0, β=0.0, γ=0.125, δ=0.125, K=[2 0 0 2])
# end
# # name: DIXMAANIGrad
# function F21(x)
#     DIXMAANHGrad(x; α=1.0, β=0.0, γ=0.125, δ=0.125, K=[2 0 0 2])
# end
# """
# Source:
# Problem 2 in Appendix

# William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
# Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493

# name: ExponetialII
# """
# function F22(x)
#     a = [i / 10 for i in 2:length(x)]
#     vcat(exp(x[1] - 1), a .* exp.(x[2:end]) .+ x[2:end] .- 1)
# end

# """
# Source:
# Problem 17 in Appendix (modified)

# William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
# Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493

# name: ExponetialIV
# """
# function F23(x)
#     n = length(x)
#     a = collect(1:length(x))
#     (a .* (exp.(x)) / n) .- 1
# end


# """
# Source:
# Problem 2 
# Wang, C., Wang, Y. & Xu, C. 
# A projection method for a system of nonlinear monotone equations with convex 
# constraints. Math Meth Oper Res 66, 33–46 (2007). 
# https://doi.org/10.1007/s00186-006-0140-y

# name: ArcTan
# """
# function F24(x)
#     n = length(x)
#     d = Normal(0.0, Float64(n))
#     M = rand(d, n, n)
#     # M = 5 * [0.33435 0.535161 0.561637 0.344291 0.111954
#     #     0.360864 0.508599 0.821533 0.694135 0.953729
#     #     0.876373 0.747615 0.336939 0.0355163 0.558138
#     #     0.983291 0.355772 0.209947 0.251682 0.410537
#     #     0.270499 0.514839 0.198478 0.883239 0.74523]
#     M = 0.5 * (M - M')
#     q = rand(d, n)
#     # q = [0.7845171302827574
#     #     0.5478931492840216
#     #     0.8419546163518192
#     #     0.3928429699008852
#     #     0.21783349212498404]
#     ρ = 100
#     q0 = 1.5 * ones(n)
#     F(y) = ρ * atan.(y .- 2) + M * y + q
#     F(x) - F(q0)
# end
# """
# Source:
# Problem 3 
# Wang, C., Wang, Y. & Xu, C. 
# A projection method for a system of nonlinear monotone equations with convex 
# constraints. Math Meth Oper Res 66, 33–46 (2007). 
# https://doi.org/10.1007/s00186-006-0140-y

# name: ArcTan2
# """
# function F25(x)
#     n = length(x)
#     d1 = Normal(0.0, 0.5)
#     d2 = Normal(0.0, 250)
#     A = rand(d1, n, n)
#     B = rand(d1, n, n)
#     B = 0.5 * (B - B')
#     M = A * A' + B
#     q = rand(d2, n)
#     q02 = 0.9 * ones(n)
#     a = rand(1:100, n)
#     D(y) = a .* atan.(y)
#     F(y) = D(y) + M * y + q
#     F(x) - F(q02)
# end

# """
# Source:
# Problem 4 
# Wang, C., Wang, Y. & Xu, C. 
# A projection method for a system of nonlinear monotone equations with convex 
# constraints. Math Meth Oper Res 66, 33–46 (2007). 
# https://doi.org/10.1007/s00186-006-0140-y

# name: StrictlyMonotone
# """
# function F26(x)
#     A = [
#         1 0 0 0
#         0 1 -1 0
#         0 1 1 0
#         0 0 0 0
#     ]
#     b = [-10; 1; -3; 0]
#     A * x + [(x[1])^3; (x[2])^3; 2(x[3])^3; 2(x[4])^3] + b
# end

# """
# Source:
# Problem 1 in 

# Abubakar, A. B., Kumam, P., & Mohammad, H. (2020). 
# A note on the spectral gradient projection method for nonlinear monotone equations with applications. 
# Computational and Applied Mathematics, 39(2), 129. https://doi.org/10.1007/s40314-020-01151-5

# name: ExponetialIIModified
# """
# function F27(x)
#     vcat(exp(x[1]) - 1, exp.(x[2:end]) .+ x[2:end] .- 1)
# end

# """
# Source:
# Problem 4.2 in 

# Liu, J., & Feng, Y. (2019). 
# A derivative-free iterative method for nonlinear monotone equations with convex constraints. 
# Numerical Algorithms, 82(1), 245–262. https://doi.org/10.1007/s11075-018-0603-2

# name: DiscreteBV
# """
# function F28(x::Vector{Float64})
#     n = length(x)
#     h = 1 / (n + 1)
#     F = zeros(n)

#     # Equation 1
#     F[1] = 2 * x[1] + 0.5 * h^2 * (x[1] + h)^3 - x[2]

#     # Equations 2 to n-1
#     for i in 2:n-1
#         F[i] = 2 * x[i] + 0.5 * h^2 * (x[i] + i * h)^3 - x[i-1] + x[i+1]
#     end

#     # Equation n
#     F[n] = 2 * x[n] + 0.5 * h^2 * (x[n] + n * h)^3 - x[n-1]

#     return F
# end

# """
# Source:
# Problem 4.3 in 

# Liu, J., & Feng, Y. (2019). 
# A derivative-free iterative method for nonlinear monotone equations with convex constraints. 
# Numerical Algorithms, 82(1), 245–262. https://doi.org/10.1007/s11075-018-0603-2

# name: TriExponential
# """
# function F29(x::Vector{Float64})
#     n = length(x)
#     F = zeros(Float64, n)

#     # Equation 1
#     F[1] = 3 * x[1]^3 + 2 * x[2] - 5 + sin(x[1] - x[2]) * sin(x[1] + x[2])

#     # Equations 2 to n-1
#     for i in 2:n-1
#         F[i] = -x[i-1] * exp(x[i-1] - x[i]) + x[i] * (4 + 3 * x[i]^2) +
#                2 * x[i+1] + sin(x[i] - x[i+1]) * sin(x[i] + x[i-1]) - 8
#     end

#     # Equation n
#     F[n] = -x[n-1] * exp(x[n-1] - x[n]) + 4 * x[n] - 3

#     return F
# end


# """
# Source:
# Problem 4.7 in 

# Liu, J., & Feng, Y. (2019). 
# A derivative-free iterative method for nonlinear monotone equations with convex constraints. 
# Numerical Algorithms, 82(1), 245–262. https://doi.org/10.1007/s11075-018-0603-2

# name: PolynomialModified
# """
# function F30(x::Vector{Float64})
#     n = length(x)
#     F = zeros(Float64, n)

#     # Equation 1
#     F[1] = 2.5 * x[1] + x[2] - 1

#     # Equations 2 to n-1
#     for i in 2:n-1
#         F[i] = x[i-1] + 2.5 * x[i] + x[i+1] - 1
#     end

#     # Equation n
#     F[n] = x[n-1] + 2.5 * x[n] - 1

#     return F
# end

# """
# Source:
# Problem 1 in 
# Cruz, W. L. (2017). 
# A spectral algorithm for large-scale systems of nonlinear monotone equations. 
# Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

# name: MinimizeMaximizeFunction
# """
# function F31(x::Vector{Float64})
#     n = length(x)
#     F = zeros(Float64, n)

#     for i in 1:n
#         F[i] = min(min(abs(x[i]), x[i]^2), max(abs(x[i]), x[i]^3))
#     end

#     return F
# end

# """
# Source:
# Problem 4 in 
# Cruz, W. L. (2017). 
# A spectral algorithm for large-scale systems of nonlinear monotone equations. 
# Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

# name: SinusoidalFunction
# """
# function F32(x::Vector{Float64})
#     n = length(x)
#     F = zeros(Float64, n)

#     F[1] = 2 * x[1] + sin(x[1]) - 1

#     for i in 2:n-1
#         F[i] = -2 * x[i-1] + 2 * x[i] + sin(x[i]) - 1
#     end

#     F[n] = 2 * x[n] + sin(x[n]) - 1

#     return F
# end



# """
# Source:
# Problem 7 in 
# Cruz, W. L. (2017). 
# A spectral algorithm for large-scale systems of nonlinear monotone equations. 
# Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

# name: LinearSystemMatrix
# """
# function F33(x::Vector{Float64})
#     n = length(x)
#     A = (5 / 2) * I(n) + diagm(-1 => ones(n - 1), 1 => ones(n - 1))
#     b = -ones(Float64, n)

#     return A * x + b
# end

# """
# Source:
# Problem 8 in 
# Cruz, W. L. (2017). 
# A spectral algorithm for large-scale systems of nonlinear monotone equations. 
# Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

# name: CubicPolynomial
# """
# function F34(x::Vector{Float64})
#     n = length(x)
#     F = zeros(Float64, n)

#     F[1] = (1 / 3) * x[1]^3 + (1 / 2) * x[2]^2

#     for i in 2:n-1
#         F[i] = -(1 / 2) * x[i]^2 + (i / 3) * x[i]^3 + (1 / 2) * x[i+1]^2
#     end

#     F[n] = -(1 / 2) * x[n]^2 + (n / 3) * x[n]^3

#     return F
# end

# """
# Source:
# Problem 9 in 
# Cruz, W. L. (2017). 
# A spectral algorithm for large-scale systems of nonlinear monotone equations. 
# Numerical Algorithms, 76(4), 1109–1130. https://doi.org/10.1007/s11075-017-0299-8

# name: QuadraticSumFunction
# """
# function F35(x::Vector{Float64})
#     n = length(x)
#     F = zeros(Float64, n)
#     smx = sum(x)
#     for i in 1:n
#         F[i] = x[i] + ((smx - x[i]^2) / n) + i
#     end

#     return F
# end




#Example 1, Γ=[-2, ∞]
function F1(x::Vector{Float64})
    n = length(x)
    result = similar(x)


    for i in 1:n
        result[i] = 2 * x[i] - sin(x[i])
    end

    return result
end
function F1v(x::Vector{Float64})
    return 2x - sin.(x)
end

#Example 2 Γ=[-1, +∞]
function F2(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = log((x[i]) + 1) - x[i] / n
    end

    return result
end
#Example 3, Γ=R^+
function F3(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = exp(x[i]) - 1
    end

    return result
end
function F4(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 1:n-1
        result[i] = 4 * x[i] + (x[i+1] - 2 * x[i]) - (x[i+1]^2) / 3
    end

    # Handle the last element separately
    result[n] = 4 * x[n] + (x[n-1] - 2 * x[n]) - (x[n-1]^2) / 3

    return result
end

#Example 6, Γ=R
function F5(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    result[1] = x[1] - exp(cos((x[1] + x[2]) / (n + 1)))
    for i in 2:n-1
        result[i] = x[i] - exp(cos((x[i-1] + x[i+1]) / (n + 1)))
    end
    result[n] = x[n] - exp(cos((x[n-1] + x[n]) / (n + 1)))

    return result
end
#Example F7, Γ=[-3, ∞]
function F6(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    result[1] = x[1] + sin(x[1]) - 1
    for i in 2:n-1
        result[i] = -x[i-1] + 2 * x[i] + sin(x[i]) - 1
    end
    result[n] = x[n] + sin(x[n]) - 1

    return result
end

function F7(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 2:n, j in 1:n
        result[1] = sum(x[j]^2)
        result[i] = -2 * x[1] * x[i]
    end

    return result
end



#Example 5

function F8(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    # Handle the first element separately
    result[1] = x[1] * (x[1]^2 + x[2]^2) - 1

    # Loop through elements 2 to n-1
    for i in 2:n-1
        result[i] = x[i] * (x[i-1]^2 + 2 * x[i]^2 + x[i+1]^2) - 1
    end

    # Handle the last element separately
    result[n] = x[n] * (x[n-1]^2 + x[n]^2)

    return result
end

function F9(x::AbstractVector{<:AbstractFloat})
    n = length(x)
    result = similar(x)


    result[1] = 3 * x[1]^3 + 2 * x[2] - 5 + sin(abs(x[1] - x[2])) * sin(abs(x[1] + x[2]))


    for i in 2:n-1
        result[i] = -x[i-1] * exp(x[i-1] - x[i]) + x[i] * (4 + 3 * x[i]^3) + 2 * x[i+1] + sin(abs(x[i] - x[i+1])) * sin(abs(x[i] + x[i+1])) - 8
    end


    result[n] = -x[n-1] * exp(x[n-1] - x[n]) + 4 * x[n] - 3

    return result
end

function F10(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = (x[i] - 1)^2 - 1.01
    end

    return result
end
function F11(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = i / n * exp(x[i]) - 1
    end

    return result
end
function F12(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = x[i] - sin(abs(x[i] - 1))
    end

    return result
end
function F13(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = 2 * x[i] - sin(abs(x[i] - 1))
    end

    return result
end
function F14(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = x[i] - 2 * sin(abs(x[i] - 1))
    end

    return result
end
function F15(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = (exp(x[i]))^2 + 3 * sin(x[i]) * cos(x[i]) - 1
    end

    return result
end


function F16(x::Vector{Float64})
    n = length(x)
    result = similar(x)

    # Handle the first element separately
    result[1] = 2.5 * x[1] + x[2] - 1

    # Loop through elements 2 to n-1
    for i in 2:n-1
        result[i] = x[i-1] + 2.5 * x[i] + x[i+1] - 1
    end

    # Handle the last element separately
    result[n] = x[n-1] + 2.5 * x[n] - 1
    return result
end
function F17(x::AbstractVector{<:AbstractFloat})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = 2 * x[i] + sin(abs(x[i]))
    end

    return result
end
function F18(x::AbstractVector{<:AbstractFloat})
    n = length(x)
    result = similar(x)

    for i in 1:n
        result[i] = 0.5 * (log(x[i]) + exp(x[i]) - sqrt((log(x[i]) - exp(x[i]))^2 - 10^(-10)))
    end

    return result
end
function F19(x::AbstractVector{<:AbstractFloat})
    n = length(x)
    result = similar(x)
    A = 0
    c = 10^(-5)
    for j in 1:n
        A += x[j]^2

    end

    for i in 1:n
        result[i] = 2 * c * (x[i] - 1) + 4 * x[i] * A - x[i]
    end

    return result
end


function F20(x::Vector{Float64})
    n = length(x)
    result = similar(x)
    A = 0
    for j in 1:n
        A += x[j]
    end

    for i in 1:n
        result[i] = x[i] * cos(x[i] - 1 / n) * (sin(x[i]) - 1 - (1 - x[i])^2 - 1 / n * A)
    end
    return result
end

"""
Li, Q., & Zheng, B. (2021). 
Scaled three-term derivative-free methods for solving large-scale nonlinear monotone equations. 
Numerical Algorithms, 87(3), 1343–1367. https://doi.org/10.1007/s11075-020-01010-8

Problem 4.8
"""

function F21(x::Vector{Float64})
    n = length(x)
    upper_diag = repeat([-1.0], n - 1)
    A = diagm(
        0 => repeat([2.0], n),
        -1 => upper_diag,
        1 => upper_diag
    )
    f = exp.(x) .- 1
    A * x + f
end

function getAllFunctions()
    return [F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21]
end

function getFunctionsWithPoints(init_points; point_labels=[:x01, :x02, :x03, :x04, :x05, :x06, :x07, :x08, :x09])
    problems = getAllFunctions()
    selected_points = map(x -> init_points[x], point_labels)
    map(x -> [x, selected_points], problems)
end