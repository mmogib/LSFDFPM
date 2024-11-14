% Example 4
function [F] = F04(x)
    n = length(x);
    F = zeros(size(x));
    for i = 1:n-1
        F(i) = 4 * x(i) + (x(i+1) - 2 * x(i)) - (x(i+1)^2) / 3;
    end
    F(n) = 4 * x(n) + (x(n-1) - 2 * x(n)) - (x(n-1)^2) / 3;
end
