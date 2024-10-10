


#Example 1, Γ=[-2, ∞]
function F1(x::Vector{Float64})
    n = length(x)
    result = similar(x)


    for i in 1:n
        result[i] = 2 * x[i] - sin(x[i])
    end

    return result
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

function getAllFunctions()
    return [F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20]
end