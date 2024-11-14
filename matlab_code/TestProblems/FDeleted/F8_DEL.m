function result = F8_DEL(x)
    n = length(x);
    result = zeros(size(x));
    result(1) = x(1) * (x(1)^2 + x(2)^2) - 1;
    for i = 2:n-1
        result(i) = x(i) * (x(i-1)^2 + 2 * x(i)^2 + x(i+1)^2) - 1;
    end
    result(n) = x(n) * (x(n-1)^2 + x(n)^2);
end