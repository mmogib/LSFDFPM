abstract type AbstractLineSearch end

struct LineSearch <: AbstractLineSearch end
struct NoLineSearch <: AbstractLineSearch end

struct LSParams
    k::Int
    Δ::Float64
    sk::Float64
    αstab::Float64
end
struct LSOutput
    α::Union{Nothing,Float64}
    iters::Int
    Δ::Float64
    αstab::Float64
end

