function result = F9_DEL(x)
    n = length(x);
    result = zeros(size(x));
    result(1) = 3 * x(1)^3 + 2 * x(2) - 5 + sin(abs(x(1) - x(2))) * sin(abs(x(1) + x(2)));
    for i = 2:n-1
        result(i) = -x(i-1) * exp(x(i-1) - x(i)) + x(i) * (4 + 3 * x(i)^3) + 2 * x(i+1) + ...
                    sin(abs(x(i) - x(i+1))) * sin(abs(x(i) + x(i+1))) - 8;
    end
    result(n) = -x(n-1) * exp(x(n-1) - x(n)) + 4 * x(n) - 3;
end
