% Example 16
function result = F16_DEL(x)
    n = length(x);
    result = zeros(size(x));
    result(1) = 2.5 * x(1) + x(2) - 1;
    for i = 2:n-1
        result(i) = x(i-1) + 2.5 * x(i) + x(i+1) - 1;
    end
    result(n) = x(n-1) + 2.5 * x(n) - 1;
end