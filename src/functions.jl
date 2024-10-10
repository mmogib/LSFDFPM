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
    xs = createInitialPoints(dim)
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