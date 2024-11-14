include("depts.jl")
include("includes.jl")


function order_line_search(r)
    O = Dict(
        "NoLineSearch" => 1,
        "LSI" => 2,
        "LSII" => 3,
        "LSIII" => 4,
        "LSIV" => 5,
        "LSV" => 6,
        "LSVI" => 7,
        "LIX" => 10,
        "LSVII" => 8,
        "LSVIII" => 9
    )
    O[r]
end

function runit(fname::String, dims::Vector{Int}, linesearchs)

    dfs_line_search = run_experiments(
        dims,
        LineSearch();
        problems_filter=i -> i <= 3,
        maxitrs=2_000,
        max_time=nothing,
        μ0=0.1,
        γ=0.2,
        ε=1e-6,
        α=0.7,
        linesearchs=linesearchs,
    )
    dfs_no_line_search = run_experiments(
        dims,
        NoLineSearch();
        problems_filter=i -> i <= 3,
        maxitrs=2_000,
        max_time=nothing,
        μ0=0.1,
        γ=0.2,
        ε=1e-6,
        α=0.7,
    )
    #     gamma = 0.2; alpha = 0.7; r=0.1; psi=0.2; gn =1; %d2
    # mu=0.1;
    DF = vcat(dfs_no_line_search, dfs_line_search)
    DF = transform(DF, :line_search => ByRow(order_line_search) => :line_search_order)
    DF = transform(DF, :problem => ByRow(x -> length(x) == 2 ? "$(x[1])0$(x[2])" : x) => :problem)
    DF = sort(DF, order(:line_search_order))
    try
        #save_plots(dfs; algorithms=Symbol.(getAllLineSearch()), label=:line_search, folder_name="./results/profiles")
        savedf("./results/dfs/$(fname)", DF, dated=true)
        save_plots(DF; algorithms=vcat(:NoLineSearch, Symbol.(linesearchs)), label=:line_search, folder_name="./results/profiles")
    catch e
        printstyled("Error $r", color=:red, bold=true)
    end
    DF
end
# 1000;  15000; 75_000; 150_000
dfs = runit("results_1k_5k_75_.xlsx", [500, 15_000, 75_000], [LSI])
# disp(dfs)
A = [13 6
    13 6
    11 5
    11 5
    10 4
    9 15
    11 16
    11 16
    11 16
    14 6
    14 6
    12 6
    12 6
    11 5
    9 15
    12 20
    12 20
    12 20
    15 7
    15 7
    13 6
    13 6
    11 5
    9 15
    13 21
    13 21
    13 21
    23 NaN
    23 NaN
    10 7
    1 NaN
    10 5
    9 28
    10 33
    10 33
    10 33
    27 NaN
    27 NaN
    11 8
    1 NaN
    11 6
    9 28
    11 37
    11 37
    11 37
    29 NaN
    28 NaN
    11 8
    1 NaN
    11 6
    9 28
    12 39
    12 39
    12 39
    3 8
    2 9
    11 6
    10 6
    10 4
    9 28
    11 33
    11 33
    11 33
    3 9
    2 9
    12 7
    11 7
    11 5
    9 28
    12 37
    12 37
    12 37
    3 9
    2 9
    13 7
    12 7
    11 5
    9 28
    12 38
    12 38
    12 38
]
df = DataFrame(
    line_search=vcat(repeat(["NoLS"], 81), repeat(["LS1"], 81)),
    iterations=vcat(A[:, 1], A[:, 2]),
    time=vcat(A[:, 1], A[:, 2]),
    func_evals=vcat(A[:, 1], A[:, 2]),
    flag="converge"
)


save_plots(df; algorithms=vcat(:NoLS, :LS1), label=:line_search, folder_name="./results/profiles")

df = XLSX.readtable("./results/dfs/very_simple.xlsx", "Sheet1") |> DataFrame


df_converged = filter(x -> occursin("converge", x[:flag]), df)
df_size = map([500, 5_000, 15_000]) do sz
    filter(x -> x[:dim] == sz, df_converged)
end
df_size[3]

p = scatter(df_size[1][!, :problem], df_size[1][!, :iterations])

# NLS_DFS_5K = XLSX.readtable("./results/dfs/with_no_line_search_s_5k_2024_10_21_08_18.xlsx", "Sheet1") |> DataFrame
# NLS_DFS_15_30K = XLSX.readtable("./results/dfs/with_no_line_search_s_15_30k_2024_10_21_09_52.xlsx", "Sheet1") |> DataFrame


# LS_DFS_5K = XLSX.readtable("./results/dfs/with_line_search_s_5k_2024_10_20_09_07.xlsx", "Sheet1") |> DataFrame
# LS_DFS_15K = XLSX.readtable("./results/dfs/with_line_search_s_15k_2024_10_20_19_26.xlsx", "Sheet1") |> DataFrame
# LS_DFS_30K = XLSX.readtable("./results/dfs/with_line_search_s_30k_2024_10_21_03_05.xlsx", "Sheet1") |> DataFrame

# DF = vcat(NLS_DFS_5K, NLS_DFS_15_30K, LS_DFS_5K, LS_DFS_15K, LS_DFS_30K)
# DF = transform(DF, :line_search => ByRow(x -> x == "compute_alpha_k" ? "LIX" : x) => :line_search)

# DF = transform(DF, :line_search => ByRow(order_line_search) => :line_search_order)
# DF = transform(DF, :problem => ByRow(x -> length(x) == 2 ? "$(x[1])0$(x[2])" : x) => :problem)
# DF = sort(DF, order(:line_search_order))

# save_plots(DF; algorithms=vcat(:NoLineSearch, Symbol.(getAllLineSearch())), label=:line_search, folder_name="./results/profiles")

# savedf("./results/dfs/ALL_DATA.xlsx", DF, dated=false)

DATA_FROM_MATLAB = XLSX.readdata("./matlab/code-withoutlinesearch/withoutLS/results/DFPMWLS_NEW11-10-2024-13-39.xlsx", "Sheet1!A1:T757")
data_header = DATA_FROM_MATLAB[1, :]
M_Iters = map(x -> ismissing(x) ? NaN : Int(x), DATA_FROM_MATLAB[2:end, 1:4:end])
MATLAB_DF = DataFrame(DF_FROM_MATLAB[2:end, :], names=DF_FROM_MATLAB[1, :])