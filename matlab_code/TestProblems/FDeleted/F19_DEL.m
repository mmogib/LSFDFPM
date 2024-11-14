% Example 19
function result = F19_DEL(x)
    n = length(x);
    A = sum(x.^2);
    c = 1e-5;
    result = 2 * c * (x - 1) + 4 * x * A - x;
end