% Example 11
function result = F08(x)
    n = length(x);
    result = (1:n)' / n .* exp(x) - 1;
end