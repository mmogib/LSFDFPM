function basenameparts(filename::String)
    parts = split(basename(filename), ".")
    if length(parts) == 1
        return parts[1], "UNKOWN"
    else
        return join(parts[1:end-1], "."), parts[end]
    end
end


function outputfilename(file_name::String; dated::Bool=true)
    fname = basename(file_name)
    root = dirname(file_name)
    root_dir = mkpath(root)
    filename = if dated
        fparts = basenameparts(fname)
        zdt = now(tz"Asia/Riyadh")
        dayfolder = Dates.format(zdt, "yyyy_mm_dd")
        hourfolder = Dates.format(zdt, "HH_MM")
        "$root_dir/$(fparts[1])_$(dayfolder)_$(hourfolder).$(fparts[2])"
    else
        "$root_dir/$fname"
    end

    filename
end

function savedf(filename::String, df::DataFrame; dated::Bool=true)
    fname = outputfilename(filename; dated)
    fparts = basenameparts(fname)
    ext = isnothing(fparts[2]) ? "UNKOWN" : lowercase(fparts[2])
    if ext in ["txt", "csv"]
        @info "Saving data in $(fname)."
        CSV.write(fname, df)
    elseif ext == "xlsx"
        @info "Saving data in $(fname)."
        XLSX.writetable(fname, df; overwrite=true)
    else
        @warn "The file $fname is not supported."
    end
end


function run_algorithm(algo::Function, dim, fns, args...; printit::Bool=false, kwargs...)
    xs = getInitialPoints(dim)
    num_rows = length(xs) * length(fns)
    T = Matrix{Any}(undef, num_rows, 8)
    counter = 0
    for f in fns
        for x_label in sort(collect(keys(xs)))
            x0_name = xs[x_label][1]
            x0 = xs[x_label][2]
            result_mopcgm = @timed algo(f, x0, args...; kwargs...)
            val = result_mopcgm.value
            t = result_mopcgm.time
            counter += 1
            T[counter, :] = [Symbol(f), x0_name, length(x0), val[1], val[2], val[3], t, val[end]]
        end
    end
    if printit
        pretty_table(T; header=["Function", "x0", "dim", "||F(x*)||",
            "Iterations", "Function Evaluations", "Time", "Flag"])
    end
    T
    df = DataFrame(
        algorithm=repeat([String(Symbol(algo))], num_rows),
        problem_name=String.(T[:, 1]),
        x0_name=String.(T[:, 2]),
        dim=T[:, 3],
        fun_norm=T[:, 4],
        iters=T[:, 5],
        fun_evaluations=T[:, 6],
        time=T[:, 7],
        message=String.(T[:, 8])
    )
    df
end

function dataFrames2Matrix(data::Vector{DataFrame}, field::Symbol, np::Int64, ns::Int64)

    T = map(data) do df
        values = df[!, field]
        values .|> x -> isempty(x) ? NaN : Float64(x)
    end |> d -> vcat(d...) |> d -> reshape(d, np, ns)
    T
end

function dataFrames2Matrix(data::Vector{DataFrame}, field::Symbol)
    ns = length(data)
    np, = size(data[1])
    dataFrames2Matrix(data, field, np, ns)
end

function readPaperData(dfs::DataFrame, names::Vector{String})
    dfs = transform(dfs, :algorithm => ByRow(strip) => :algorithm)
    data = map(enumerate(names)) do (i, name)
        filter(row -> row[:algorithm] == name, dfs)
    end
    data

end
function readPaperData(input_file::String, names::Vector{String})
    dfs = CSV.File(input_file) |> DataFrame
    readPaperData(dfs, names)
end

function plotPerformanceProfile(T::Matrix{<:Number}, title::String, labels::Vector{String}; colors::Union{ColorPalette,Vector{Symbol}}=palette(:tab10), logscale::Bool=true)
    (w, h) = Plots._plot_defaults[:size]
    (_, _, max_ratio) = performance_profile_data(T)
    p = performance_profile(PlotsBackend(),
        T, labels;
        size=(1.2w, h),
        logscale,
        title=title,
        xlabel=L"\tau",
        ylabel=L"\rho(\tau)",
        legendfontsize=8,
        linestyle=:dash,
        palette=colors,
        linewidth=2.5,
        minorgrid=true, leg=:bottomright
    )
    p = if (logscale)
        xs = 1:ceil(max_ratio + 0.5)
        ys = 0:0.1:1.01
        plot(p,
            xticks=(xs, map(x -> "$(Int(x))", xs)),
            yticks=(ys, map(x -> "$x", ys)),
            # framestyle=:origin
        )
    else
        p
    end
    p

end
function produceProfileImages(input_file::String)
    colors = [:red, :black, :green, :blue, :purple]

    output_folder = "./results/profiles/imgs"


    names = String.(Symbol.([MOPCGM, GMOPCGM, CGPM, GCGPM, STTDFPM, Framework]))
    printstyled("reading stored data in "; color=:blue)
    printstyled("$input_file\n"; color=:blue, bold=true)
    data = readPaperData(input_file, names)
    printstyled("creating plots ...\n"; color=:green)

    ns = length(data)
    np, = size(data[1])
    # 
    plts = map([
        (:iters, "Iterations"),
        (:fun_evaluations, "Function Evaluations"),
        (:time, "CPU Time")]) do (item, title)
        T = dataFrames2Matrix(data, item, np, ns)
        # pretty_table(T)
        printstyled("creating plots for $title \n"; color=:green)
        p = plotPerformanceProfile(T, title, names; colors)
        p, title
    end

    map(plts) do (p, title)
        file_name = replace(title, " " => "_") |> lowercase
        printstyled("saving plots for $title to file \n"; color=:reverse)
        png_file = outputfilename("$output_folder/$(file_name).png"; dated=false)
        savefig(p, png_file)
        # svg_file = outputfilename(file_name, "svg"; root=output_folder, suffix="performance_profiles")
        # savefig(p, svg_file)
        # pdf_file = outputfilename(file_name, "pdf"; root=output_folder, suffix="performance_profiles")
        # savefig(p, pdf_file)
        # eps_file = outputfilename(file_name, "eps"; root=output_folder, suffix="performance_profiles")
        # savefig(p, eps_file)
    end
    printstyled("Finished saving in $output_folder\n"; color=:green)
    plts
end

function save_plots(df::DataFrame; algorithms::Vector{Symbol}, label::Symbol, folder_name::String)
    colors = [:red,          # 1
        :blue,         # 2
        :green,        # 3
        :darkorange1,  # 4
        :purple,       # 5
        :yellow,       # 6
        :cyan,         # 7
        :magenta,      # 8
        :black,        # 9
        :lime
    ]
    (w, h) = Plots._plot_defaults[:size]
    ys = 0:0.1:1.01

    plts = map([(:iterations, "Iterations"), (:time, "Time"),
        (:func_evals, "Function Evaluations")]) do (data_point, data_header)
        T = Matrix{Float64}(undef, Int(DataFrames.nrow(df) / length(algorithms)), length(algorithms))
        foreach(enumerate(algorithms)) do (i, name)
            d = select(filter(r -> r[label] == String(name), df),
                [data_point, :flag] => ByRow((itr, f) -> occursin("converge", f) ? itr : NaN) => Symbol(data_header))
            T[:, i] = d[!, Symbol(data_header)]
        end
        pdata = performance_profile_data(T)
        max_ratio = pdata[3]
        xs = 1:ceil(max_ratio + 0.5)
        p = performance_profile(PlotsBackend(), T, String.(algorithms);
            title="$(data_header)",
            logscale=true,
            size=(1.2w, h),
            xlabel="Performance ratio",
            ylabel="Solved problems (%)",
            legendfontsize=8,
            linestyle=:dash,
            palette=colors,
            linewidth=2.5,
            minorgrid=true, leg=:bottomright
        )
        plot(p, xticks=(xs, map(x -> "$(Int(x))", xs)),
            yticks=(ys, map(x -> "$x", ys))), data_header
    end

    for (p, name) in plts
        # output_file_svg = outputfilename("$(folder_name)_$name", "svg"; suffix)
        output_file_png = outputfilename("$(folder_name)/$(name).png")

        # Plots.svg(p, output_file_svg)
        Plots.png(p, output_file_png)
    end
end



function run_experiments(
    dims::Vector{Int},
    LS::AbstractLineSearch;
    problems_filter::Function=i -> i != 9,
    maxitrs=2_000,
    max_time::Union{Nothing,Integer}=nothing,
    μ0=0.5,
    γ=1.8,
    ε=1e-6,
    α=0.7,
    linesearchs=[LIX],
)
    search_directions = mohammed_ibrahim_suliman #getAllSearchDirection()
    lss = if isa(LS, NoLineSearch)
        nothing
    else
        linesearchs
    end

    io = open("./results/error_log.txt", "w+")
    logger = SimpleLogger(io)
    dfs = with_logger(logger) do
        return run_experiments(dims, search_directions, lss;
            problems_filter=problems_filter,
            maxitrs=maxitrs,
            max_time=max_time,
            μ0=μ0,
            γ=γ,
            ε=ε,
            α=α,
        )
    end
    close(io)
    dfs

end

function run_experiments(dims,
    search_directions,
    linesearches;
    problems_filter=i -> i != 9,
    maxitrs=2_000,
    max_time::Union{Nothing,Integer}=nothing,
    μ0=0.5,
    γ=1.8,
    ε=1e-6,
    α=0.7,
)

    log_print(sd, ls, f, x0_name, dim, sol) = isnothing(linesearches) ? join([String(Symbol(sd)), " : ",
            "NoLineSearch", " : ",
            String(Symbol(f)), ":",
            x0_name, " : ",
            "dim: ", dim, " : ",
            sol, "\n"], "") : join([String(Symbol(sd)), " : ",
            String(Symbol(ls)), " : ",
            String(Symbol(f)), ":",
            x0_name, " : ",
            "dim: ", dim, " : ",
            sol, "\n"], "")


    lss = isnothing(linesearches) ? [(g, m) -> nothing] : isa(linesearches, Vector) ? linesearches : [linesearches]
    results_sd = map([search_directions]) do sd
        df_lss = map(lss) do ls
            df_dim = map(dims) do dim
                all_points = getInitialPoints(dim)
                problems = getFunctionsWithPoints(all_points)
                problems = problems[map(problems_filter, 1:length(problems))]
                df_f = map(problems) do (f, points)
                    df_p = map(points) do (x0_label, x0_name, x0)
                        df = DataFrame(
                            search_direction=String(Symbol(sd)),
                            line_search=isnothing(linesearches) ? "NoLineSearch" : String(Symbol(ls)),
                            problem=String(Symbol(f)),
                            dim=dim,
                            x0_label=String(x0_label),
                            x0_name=x0_name,
                            norm_f=NaN,
                            iterations=NaN,
                            func_evals=NaN,
                            time=NaN,
                            flag="Error"
                        )
                        try
                            t_sol = @timed LSFDFPM(
                                f,
                                x0,
                                sd(),
                                ls(f, maxitrs),
                                maxitrs,
                                ε,
                                max_time,
                                μ0,
                                γ,
                                α
                            )

                            t = t_sol.time
                            sol = t_sol.value
                            norm_f = sol[1]
                            iterations = sol[2]
                            func_evals = sol[3]
                            flag = String(sol[end])

                            df = DataFrame(
                                search_direction=String(Symbol(sd)),
                                line_search=isnothing(linesearches) ? "NoLineSearch" : String(Symbol(ls)),
                                problem=String(Symbol(f)),
                                dim=dim,
                                x0_label=String(x0_label),
                                x0_name=String(x0_name),
                                norm_f=norm_f,
                                iterations=iterations,
                                func_evals=func_evals,
                                time=t,
                                flag=flag
                            )
                            c = occursin("converge", flag) ? :green : :blue
                            printstyled(log_print(sd, ls, f, x0_name, dim, sol), color=c)
                        catch e
                            printstyled(log_print(sd, ls, f, x0_name, dim, "Error"), color=:red)

                            timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
                            @error "[$timestamp] An error occured: $e"
                        end
                        df
                    end
                    vcat(df_p...)
                end
                vcat(df_f...)
            end
            vcat(df_dim...)
        end
        vcat(df_lss...)
    end
    results_df = vcat(results_sd...)
    results_df
end