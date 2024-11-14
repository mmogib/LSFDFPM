% Example 18
function result = F14(x)
    result = 0.5 * (log(x) + exp(x) - sqrt((log(x) - exp(x)).^2 - 1e-10));
end