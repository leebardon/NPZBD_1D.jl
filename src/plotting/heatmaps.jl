using NCDatasets
using Plots, ColorSchemes
using DataFrames
using SparseArrays, LinearAlgebra, Statistics

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")

function final_year(vars, ds)
    # 100 ts each day, 1 in 5 recorded (i.e. 20 each day)
    # final year = 366 * 20 = 7320 recorded ts

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


function plot_heatmaps(fsaven, var)

    zc = get_zc(400)
    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")

    p, b = final_year(["p", "b"], ds)

    # heatmaps(p, zc, filename, "P", var)
    heatmaps(b, zc, filename, "B", var)

    # for species in [p, b]
    #     heatmaps(species, zc, ds, filename, species_names)
    # end

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


function heatmaps(species, zc, filename, s_name, var)

    parent_folder = "results/plots/heatmaps/"
    dir = check_subfolder_exists(filename, parent_folder)

    d = -zc[1:40]
    num_species = size(species, 2)
    fig = Array{Plots.Plot, 1}(undef, num_species)

    # d1 = -zc[1:20]
    # d2 = -zc[20:40]
    # fig1 = Array{Plots.Plot, 1}(undef, num_species)
    # fig2 = Array{Plots.Plot, 1}(undef, num_species)

    for i in range(1, num_species)

        data = species[1:40, i, :]
        survivors = set_extinct_to_zero(data)
        daily_data = survivors[:, 1:20:end]
        days = collect(1:size(daily_data, 2))

        z =  get_z_axis(d, days, daily_data)
        fig[i] = heatmap(days, reverse(d), reverse(z), xrotation=45, clim=(0, 0.28))
        println("done $i")

        # data_top = collect(eachslice(species[1:20, i, :], dims=2))
        # daily_top = data_top[1:20:end]
        # data_mid = collect(eachslice(species[20:40, i, :], dims=2))
        # daily_mid = data_mid[1:20:end]
        # days = collect(1:(length(daily_top)))
        # ztop =  get_z_axis(d1, days, daily_top)
        # zmid = get_z_axis(d2, days, daily_mid)
        # fig1[i] = heatmap(days, reverse(d1), reverse(ztop), xrotation=45)
        # fig2[i] = heatmap(days, reverse(d2), reverse(zmid), xrotation=45)

    end

    f = plot(fig..., 
    fg_legend = :transparent,
    size=(1300,1000),
    plot_title = "$var over time (0-400m)"
    )

    savefig(f, "$(dir)/$(s_name)_$(filename).png")


    # f1 = plot(fig1..., 
    # fg_legend = :transparent,
    # # layout = (2,nd),
    # size=(1300,1000),
    # plot_title = "Surface (0-200m)"
    # )
    # f2 = plot(fig2..., 
    # fg_legend = :transparent,
    # # layout = (2,nd),
    # size=(1300,1000),
    # plot_title = "Mid-Depths (200-400m)"
    # )

    # savefig(f1, "$(dir)/top_$(s_name)_$(filename).png")
    # savefig(f2, "$(dir)/mid_$(s_name)_$(filename).png")

end



fsaven = "results/outfiles/Wi100y_231011_23:28_8P20Z13B5D.nc"
plot_heatmaps(fsaven, "Biomass")