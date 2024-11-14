% Example 20
function result = F20_DEL(x)
    n = length(x);
    A = sum(x);
    result = x .* cos(x - 1 / n) .* (sin(x) - 1 - (1 - x).^2 - 1 / n * A);
end