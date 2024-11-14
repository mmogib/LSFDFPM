function result = F05(x)
    n = length(x);
    result = zeros(size(x));
    result(1) = x(1) - exp(cos((x(1) + x(2)) / (n + 1)));
    for i = 2:n-1
        result(i) = x(i) - exp(cos((x(i-1) + x(i+1)) / (n + 1)));
    end
    result(n) = x(n) - exp(cos((x(n-1) + x(n)) / (n + 1)));
end