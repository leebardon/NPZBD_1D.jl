using NCDatasets
# using Plots, ColorSchemes
using CairoMakie
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


function plot_bmass_heatmaps(fsaven, varname)

    zc = get_zc(400)
    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")

    p, b = final_year(["p", "b"], ds)

    varname == "P" ? bmass_heatmaps(p, zc, filename, varname) : bmass_heatmaps(b, zc, filename, varname)

end


function get_z_axis(depth, days, daily_data)

    z = Array{Float64}(undef, size(depth, 1), size(daily_data, 2))

    for i in eachindex(depth)
        for j in eachindex(days)
            z[i, j] = daily_data[i, j]
        end
    end

    return z'

end


function set_zmax(species, num_species)

    zmax = 0
    for s in range(1, num_species)
        species_max = maximum(species[1:40, s, :])
        species_max > zmax ? zmax = species_max : nothing
    end

    return zmax

end


function bmass_heatmaps(species, zc, filename, varname)

    parent_folder = "results/plots/heatmaps/biomass/"
    dir = check_subfolder_exists(filename, parent_folder)
    
    varname == "P" ? max_d = 30 : max_d = 30
    d = -zc[1:max_d]
    num_species = size(species, 2)
    zmax = set_zmax(species, num_species)
    println(zmax)
    lfs=9
    varname == "P" ? res=(1000, 500) : res=(1500, 1200)

    joint_limits = (0, zmax)
    row = 1
    col = 1
    fig = Figure(resolution=res)
    for i in range(1, num_species)
        data = species[1:max_d, i, :]
        survivors = set_extinct_to_zero(data)
        daily_data = survivors[:, 1:20:end]
        days = collect(1:size(daily_data, 2))
        z = get_z_axis(d, days, daily_data)

        heatmap(fig[row, col], days, reverse(d), reverse(z), xrotation=45,  colorrange=joint_limits)
        col += 1
        if mod(i, 5) == 0 
            row += 1
            col -= 5
        end
    end

    Colorbar(fig[:, end+1], colorrange=joint_limits, size=30, label=L"mmol/m^3", labelsize=20)

    save("$(dir)/$(varname)_$(filename).png", fig)
end


function copio_heatmaps(fsaven, copio, s_name)

    parent_folder = "results/plots/heatmaps/copio/"
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    zc = get_zc(400)
    depth = -zc[1:30]

    data = copio[1:30, :]
    daily_data = data[:, 1:20:end]
    days = collect(1:size(daily_data, 2))
 
    z = daily_data'
    joint_limits = (minimum(daily_data), maximum(daily_data))
    fig = Figure(resolution=(600,500))
    ax, hm = heatmap(fig[1, 1], days, reverse(depth), reverse(z), clim=joint_limits, colormap=(:berlin50))

    ax.xlabel="Days"
    ax.ylabel="Depth (m)"
    ax.title="$s_name Copiotrophy Index"

    Colorbar(fig[:, end+1], colorrange=joint_limits, colormap=(:berlin50), size=20)

    println("Saving fig to $(dir)/$(s_name)_$(filename).png")
    save("$(dir)/$(s_name)_$(filename).png", fig)

end




# fsaven = "results/outfiles/Wi100y_231011_20:23_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi100y_231011_23:28_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi100y_231017_01:23_10P3Z21B9D.nc"
# fsaven = "results/outfiles/Wi100y_231017_12:01_10P3Z21B9D.nc"
# fsaven = "results/outfiles/Wi50y_231017_23:00_10P3Z18B8D.nc"
# fsaven = "results/outfiles/Wi50y_231017_23:32_10P3Z18B8D.nc"
# plot_bmass_heatmaps(fsaven, "P")
# plot_bmass_heatmaps(fsaven, "B")