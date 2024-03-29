# using NCDatasets, DataFrames
# using CairoMakie
# using LinearAlgebra, Statistics

# include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
# include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")



function plot_bmass_heatmaps(fsaven, varname)

    zc = get_zc(890)
    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")

    # p, b = get_final_year(ds, ["p", "b"])

    ds.attrib["Season"] == "winter" ? pulse_freq = 10 : pulse_freq = 30
    p, b = get_final_three_cycles(ds, ["p", "b"], pulse_freq)
    varname == "P" ? bmass_heatmaps(p, zc, filename, varname) : bmass_heatmaps(b, zc, filename, varname)

end


function bmass_heatmaps(state_var, zc, filename, varname)

    parent_folder = "results/plots/heatmaps/biomass/"
    dir = check_subfolder_exists(filename, parent_folder)

    varname == "P" ? max_d = 30 : max_d = 89
    d = -zc[1:max_d]
    num_state_var = size(state_var, 2)
    zmax, z97 = set_zmax(state_var, num_state_var)
    println(zmax)
    lfs=9
    varname == "P" ? res=(1000, 500) : res=(1200, 1800)

    joint_limits = (0, zmax)
    row, col = 1, 1
    fig = Figure(size=res)
    for i in range(1, num_state_var)
        data = state_var[1:max_d, i, :]
        survivors = set_extinct_to_zero(data)
        daily_data = survivors[:, 1:20:end]
        days = collect(1:size(daily_data, 2))
        z = get_hmap_z_axis(d, days, daily_data)

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


function npp_heatmaps(fsaven, npp_months, npp, mov_avs, daily_means)
       
    parent_folder = "results/plots/heatmaps/npp/"
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    cols_hm = :tofino25
    cols_ln = ["orangered", "seagreen3", "purple4", "black", "skyblue1"]
    dnames = ["5m", "50m", "100m", "150m", "200m"]

    zc = get_zc(890)
    depth = -zc[1:30]
    ts = collect(1:length(mov_avs[1]))

    let
        fig = Figure(resolution=(600,1250))

        z = npp'
        limits = (minimum(npp), maximum(npp))
        ax1, hm = heatmap(fig[1:2, 1], ts, reverse(depth), reverse(z), clim=limits, colormap=(cols_hm))
        Colorbar(fig[1:2, 2], colorrange=limits, colormap=(cols_hm), size=20, label=L"mmol/m^3", labelsize=20)
        ax1.xlabel=L"Timesteps"
        ax1.ylabel=L"Depth~(m)"
        ax1.title="Net Primary Productivity"

        z2 = npp_months'
        days = collect(1:size(npp_months, 2))
        limits = (minimum(npp_months), maximum(npp_months))
        ax2, hm = heatmap(fig[3:4, 1], days, reverse(depth), reverse(z2), clim=limits, colormap=(cols_hm))
        Colorbar(fig[3:4, 2], colorrange=limits, colormap=(cols_hm), size=20, label=L"mmol/m^3", labelsize=20)
        ax2.xlabel=L"Months"
        ax2.xticks = 1:12
        ax2.ylabel=L"Depth~(m)"
        ax2.title="NPP Monthly Means"

        ax3 = fig[5, 1] = Axis(fig)
        for i in range(1, length(dnames))
            lines!(ts, mov_avs[i], alpha=0.7, linewidth=3, color=cols_ln[i], label=dnames[i])
             # td = collect(1:length(daily_means[1]))
            # scatterlines!(t, convert(Array{Float64}, daily_means[i]), markersize=3, color=cols_ln[i], strokewidth=1, strokecolor=:black, alpha=0.1, label=dnames[i])
            # lines!(td, convert(Array{Float64}, daily_means[i]), color=cols_ln[i], alpha=0.5, label=dnames[i], grid=false)
        end

        fig[5, 2] = Legend(fig, ax3, framevisible=false)
        ax3.xlabel=L"Timesteps"
        ax3.ylabel=L"mmol/m^3"
        ax3.title="Moving Averages"

        println("Saving fig to $(dir)/$(filename)_npp.png")
        save("$(dir)/$(filename)_npp.png", fig)
    end

end


# fsaven = "results/outfiles/Wi50y_231122_16:59_6P3Z13B8D.nc"
# fsaven = "results/outfiles/Su50y_231122_20:39_6P3Z13B8D.nc"
# fsaven = "results/outfiles/Wi50y_231122_22:17_6P3Z13B8D.nc"

# fsaven = "results/outfiles/Su50y_240104_14:04_6P3Z13B8D.nc"
# fsaven = "results/outfiles/Wi50y_240103_14:18_6P3Z13B8D.nc"

# fsaven="results/outfiles/Wi50y_240108_16:48_6P3Z13B8D.nc"
# fsaven="results/outfiles/Wi50y_240108_20:47_6P3Z13B8D.nc"
# plot_bmass_heatmaps(fsaven, "P")
# plot_bmass_heatmaps(fsaven, "B")


# plot_bmass_heatmaps(fsaven, "B")

# for var in ["P", "B"]
#     plot_bmass_heatmaps(fsaven, var, bloom)
# end