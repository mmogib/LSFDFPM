function x0 = getInitialPoint(xnum, dim)
    switch xnum
    case 1
        x0 = 10 * ones(dim, 1);  % All elements are 10
    case 2
        x0 = -10 * ones(dim, 1);  % All elements are -10
    case 3
        x0 = ones(dim, 1);  % All elements are 1
    case 4
        x0 = -ones(dim, 1);  % All elements are -1
    case 5
        x0 = 0.1 * ones(dim, 1);  % All elements are 0.1
    case 6
        x0 = 1 ./ (1:dim)';  % Elements are 1/i (vectorized)
    case 7
        x0 = (1:dim)' / dim;  % Elements are i/dim (vectorized)
    case 8
        x0 = 1 - (1:dim)' / dim;  % Elements are 1 - i/dim (vectorized)
    case 9
        x0 = rand(dim, 1);  % Random elements in [0, 1]
        otherwise, x0 = xnum;  % Use the provided value directly
    end
end
