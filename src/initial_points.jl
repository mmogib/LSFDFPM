
# Initial Point x0
# x1 = 0.0 * ones(500);
# x2 = 0.2 * ones(500);
# x3 = 0.4 * ones(500);
# x4 = 0.6 * ones(500);
# x5 = 0.8 * ones(500);
# x6 = 1.0 * ones(500);
# x7 = 1.1 * ones(500);
# x8 = 8.0 * ones(500);

function createInitialPoints(dim::Integer)
    onez = ones(dim)
    indxs = 1:dim
    startin_points_1 = [
        (:x01, "0.0", zeros(dim)),
        (:x02, "0.2", 8.2 * onez),
        (:x03, "0.4", 0.4 * onez),
        (:x04, "0.5", 0.5 * onez),
        (:x05, "0.6", 0.6 * onez),
        (:x06, "0.8", 0.8 * onez),
        (:x07, "1.0", onez),
        (:x08, "1.1", 1.1 * onez),
    ]
    startin_points_2 = [
       (:x09, "1 - 1/n", indxs .|> t -> 1 - (1 / dim)),
       (:x10, "1/k", indxs .|> t -> 1 / t),
       (:x11, "(k-1)/n", indxs .|> t -> (t - 1) / dim),
       (:x12, "1/n", indxs .|> t -> 1 / dim),
       (:x13, "1/3ᵏ", indxs .|> t -> 1 / BigFloat(3^t) .|> Float64),
       (:x14, "k/n", indxs .|> t -> t / dim)
    ]

    startin_points = [
        startin_points_1..., startin_points_2...
    ]
    # startin_points = [
    #     ("x1: 1", onez),
    #     ("x2: 0.2", 0.2 * onez),
    #     ("x3: 1/2ᵏ", indxs .|> t -> 1 / 2^t),
    #     ("x4: (k-1)/n", indxs .|> t -> (t - 1) / dim),
    #     ("x5: 1/k", indxs .|> t -> 1 / t),
    #     ("x6: 1/n", indxs .|> t -> 1 / dim),
    #     ("x7: 1 - 1/n", indxs .|> t -> 1 - (1 / dim)),
    #     ("x8: 1.1", 1.1 * onez),
    # ]
    return Dict(map(x -> x[1] => x[2:3], startin_points))
end
