function result = F06(x)
    n = length(x);
    result = zeros(size(x));
    result(1) = x(1) + sin(x(1)) - 1;
    for i = 2:n-1
        result(i) = -x(i-1) + 2 * x(i) + sin(x(i)) - 1;
    end
    result(n) = x(n) + sin(x(n)) - 1;
end
