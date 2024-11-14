
# Initial Point x0
# x1 = 0.0 * ones(500);
# x2 = 0.2 * ones(500);
# x3 = 0.4 * ones(500);
# x4 = 0.6 * ones(500);
# x5 = 0.8 * ones(500);
# x6 = 1.0 * ones(500);
# x7 = 1.1 * ones(500);
# x8 = 8.0 * ones(500);

function getInitialPoints(dim::Integer)
    onez = ones(dim)
    indxs = 1:dim
    # case 1
    #     x0 = 10 * ones(dim, 1);  % All elements are 10
    # case 2
    #     x0 = -10 * ones(dim, 1);  % All elements are -10
    # case 3
    #     x0 = ones(dim, 1);  % All elements are 1
    # case 4
    #     x0 = -ones(dim, 1);  % All elements are -1
    # case 5
    #     x0 = 0.1 * ones(dim, 1);  % All elements are 0.1
    # case 6
    #     x0 = 1 ./ (1:dim)';  % Elements are 1/i (vectorized)
    # case 7
    #     x0 = (1:dim)' / dim;  % Elements are i/dim (vectorized)
    # case 8
    #     x0 = 1 - (1:dim)' / dim;  % Elements are 1 - i/dim (vectorized)
    # case 9
    #     x0 = rand(dim, 1);  % Random elements in [0, 1]
    Random.seed!(2024)
    N_d = truncated(Normal(), 0.1, 0.9999)
    startin_points = [
        (:x01, "10", 10.0 * onez),
        (:x02, "-10", -10.0 * onez),
        (:x03, "1.0", onez),
        (:x04, "-1.0", -onez),
        (:x05, "0.1", 0.1 * onez),
        (:x06, "1/k", indxs .|> t -> 1 / t),
        (:x07, "k/n", indxs .|> t -> t / dim),
        (:x08, "(1-k)/n", indxs .|> t -> (1 - t) / dim),
        (:x09, "random", rand(N_d, dim)),
        (:x10, "1/n", indxs .|> t -> 1 / dim),
        (:x11, "0.0", zeros(dim)),
        (:x12, "0.4", 0.4 * onez),
        (:x13, "0.5", 0.5 * onez),
        (:x14, "0.6", 0.6 * onez),
        (:x15, "0.8", 0.8 * onez),
        (:x16, "1.1", 1.1 * onez),
        (:x17, "5.0", 5.0 * onez),
        (:x18, "1 - 1/n", indxs .|> t -> 1 - (1 / dim)),
        (:x19, "(k-1)/n", indxs .|> t -> (t - 1) / dim),
        (:x20, "1/3ᵏ", indxs .|> t -> 1 / BigFloat(3^t) .|> Float64),
    ]

    return Dict(map(x -> x[1] => x, startin_points))
end
