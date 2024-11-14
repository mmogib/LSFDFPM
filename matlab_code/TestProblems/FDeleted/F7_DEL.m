function result = F7_DEL(x)
    n = length(x);
    result = zeros(size(x));
    result(1) = sum(x.^2);
    for i = 2:n
        result(i) = -2 * x(1) * x(i);
    end
end