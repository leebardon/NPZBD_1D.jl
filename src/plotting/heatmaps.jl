using NCDatasets
using Plots, ColorSchemes
using DataFrames
using SparseArrays, LinearAlgebra, Statistics

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")


function final_year(vars, ds)
    # 100 ts each day, 1 in 5 recorded (i.e. 20 each day) -- 366 * 20 = 7320 recorded ts

    final_yr = Vector{Any}()

    for v in vars
        if v != "o"
            append!(final_yr, [ds[v][:, :, end-7319:end]])
        else
            append!(final_yr, [ds[v][:, end-7319:end]])
        end
    end

    return final_yr

end


function plot_bmass_heatmaps(fsaven, var)

    zc = get_zc(400)
    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")

    p, b = final_year(["p", "b"], ds)

    bmass_heatmaps(p, zc, filename, "P", var)
    # heatmaps(b, zc, filename, "B", var)

end


function get_z_axis(depth, days, daily_data)

    z = Array{Float64}(undef, size(depth, 1), size(daily_data, 2))

    for i in eachindex(depth)
        for j in eachindex(days)
            z[i, j] = daily_data[i, j]
        end
    end

    return z

end


function set_zmax(species, num_species)

    zmax = 0
    for s in range(1, num_species)
        species_max = maximum(species[1:40, s, :])
        species_max > zmax ? zmax = species_max : nothing
    end

    return zmax

end


function bmass_heatmaps(species, zc, filename, s_name, var)

    parent_folder = "results/plots/heatmaps/biomass/"
    dir = check_subfolder_exists(filename, parent_folder)
    
    d = -zc[1:40]
    num_species = size(species, 2)
    zmax = set_zmax(species, num_species)

    fig = Array{Plots.Plot, 1}(undef, num_species)

    for i in range(1, num_species)
        data = species[1:40, i, :]
        survivors = set_extinct_to_zero(data)
        daily_data = survivors[:, 1:20:end]
        days = collect(1:size(daily_data, 2))
        z = get_z_axis(d, days, daily_data)
        fig[i] = heatmap(days, reverse(d), reverse(z), xrotation=45, clim=(0, zmax))
    end

    f = plot(fig..., 
    fg_legend = :transparent,
    size=(1300,1000),
    plot_title = "$s_name $var over time (0-400m)"
    )

    println("Saving fig to $(dir)/$(s_name)_$(filename).png")
    savefig(f, "$(dir)/$(s_name)_$(filename).png")

end


function copio_heatmaps(fsaven, copio, s_name)

    parent_folder = "results/plots/heatmaps/copio/"
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    zc = get_zc(400)

    if s_name == "P" 
        depth = -zc[1:25] 
        data = copio[1:25, :]
    else
        depth = -zc
        data = copio[1:40, :]
    end 

    daily_data = data[:, 1:20:end]
    days = collect(1:size(daily_data, 2))
 
    fig = heatmap(days, reverse(depth), reverse(daily_data), xrotation=45, clim=(0.3, 0.8), c=:berlin50,
    xlabel="Days", ylabel="Depth (m)", title="$s_name Copiotrophy Index")

    println("Saving fig to $(dir)/$(s_name)_$(filename).png")
    savefig(fig, "$(dir)/$(s_name)_$(filename).png")

end




# fsaven = "results/outfiles/Wi100y_231011_20:23_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi100y_231011_23:28_8P20Z13B5D.nc"
# plot_heatmaps(fsaven, "biomass")