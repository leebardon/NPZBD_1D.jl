using CSV, Plots, DataFrames, LaTeXStrings, StatsPlots
using Plots.PlotMeasures, Statistics

include("plotting/spot_data_plots.jl")



function plot_bottle_data(colname)

    btl_means, bwi, bsp, bsu, bfa = get_btl_data()

    wi_means = filter(row -> row.season == "winter", btl_means)
    sp_means = filter(row -> row.season == "spring", btl_means)
    su_means = filter(row -> row.season == "summer", btl_means)
    fa_means = filter(row -> row.season == "fall", btl_means)
    overall_mean = combine(groupby(all, :CTDDEPTH), colname => mean => :all_mean)

    # for ds in [wi, sp, su, fa]
    #     replace!(ds[!,colname], NaN=>missing)
    #     dropmissing!(ds, :CTDTMP)
    # end

    # wi_means = filter(row -> row.season == "winter", btl_means)
    # sp_means = filter(row -> row.season == "spring", btl_means)
    # su_means = filter(row -> row.season == "summer", btl_means)
    # fa_means = filter(row -> row.season == "fall", btl_means)

    # p1 = scatter(win[:,5], -win[:,4], markersize=7, alpha=0.7, label=false, grid=false, border=:box, title="SPOT CDT Temperature", 
    # ylabel="Depth(m)", xlabel=xlab)
    # scatter!(spr[:,5], -spr[:,4], markersize=7, alpha=0.7, label=false)
    # scatter!(su[:,5], -su[:,4], markersize=7, alpha=0.7, label=false)
    # scatter!(fa[:,5], -fa[:,4], markersize=7, alpha=0.7, label=false)
    # boxplot!(x, y, grid=false, label=false, line=(1, :black), alpha=0.6,inset=(1, bbox(0.05, 0.2, 0.5, 0.25, :bottom, :right)), 
    # subplot=2, title="DCM Depth", border=:box)
    # savefig(p1, "dcm.png")

end


function btl_boxplots()

    DCM_win = filter(row -> row.DepthBin == "DCM", win)
    DCM_spr = filter(row -> row.DepthBin == "DCM", spr)
    DCM_sum = filter(row -> row.DepthBin == "DCM", sum)
    DCM_fall = filter(row -> row.DepthBin == "DCM", fall)
    # remove missing
    DCM_wi = dropmissing(DCM_win, :depth)
    DCM_sp = dropmissing(DCM_spr, :depth)
    DCM_su = dropmissing(DCM_sum, :depth)
    DCM_fa = dropmissing(DCM_fall, :depth)

    x = [["Winter"], ["Spring"], ["Summer"], ["Fall"]]
    y = [[-DCM_wi[:,4]], [-DCM_sp[:,4]], [-DCM_su[:,4]], [-DCM_fa[:,4]]]

    boxplot!(x[1], y[1], grid=false, label=false, line=(1, :black), fill=(0.3, :blue), inset=(1, bbox(0.05, 0.1, 0.5, 0.25, :bottom, :right)), subplot=2)
    boxplot!(x[2], y[2], grid=false, label=false, line=(1, :black), fill=(0.3, :green), inset=(1, bbox(0.05, 0.1, 0.5, 0.25, :bottom, :right)), subplot=2)
    boxplot!(x[3], y[3], grid=false, label=false, line=(1, :black), fill=(0.3, :red), inset=(1, bbox(0.05, 0.1, 0.5, 0.25, :bottom, :right)), subplot=2)
    boxplot!(x[4], y[4], grid=false, label=false, line=(1, :black), fill=(0.3, :orange), inset=(1, bbox(0.05, 0.1, 0.5, 0.25, :bottom, :right)), subplot=2)

end


function plot_ctd_data(maxdepth=0, varname="Temperature", colname="CTDTMP", xlab=L"\degree C")

    ctd_means, wi, sp, su, fa, all = get_ctd_data(maxdepth)
    wi_means = filter(row -> row.season == "winter", ctd_means)
    sp_means = filter(row -> row.season == "spring", ctd_means)
    su_means = filter(row -> row.season == "summer", ctd_means)
    fa_means = filter(row -> row.season == "fall", ctd_means)
    overall_mean = combine(groupby(all, :CTDDEPTH), colname => mean => :all_mean)

    for ds in [wi, sp, su, fa]
        replace!(ds[!,colname], NaN=>missing)
        dropmissing!(ds, colname)
    end

    p1, p2, p3, p4 = plot_hist2d(wi, sp, su, fa, overall_mean, xlab, colname)
    p5 = ridgeline_plot([wi[!,colname], sp[!,colname], su[!,colname], fa[!,colname]], xlab, colname)

    if maxdepth == 0
        p6, _ = plot_line_box_violin(wi, sp, su, fa, wi_means, sp_means, su_means, fa_means, xlab, colname, maxdepth)
        l = @layout [ a{0.5w} [grid(2, 2)] 
                            b{0.5h}]
                            
        fig = plot(p5, p1, p2, p3, p4, p6,
            layout = l,
            fg_legend = :transparent,
            size=(900,900),
            plot_title = "Seasonal CTD $(varname) at SPOT",
        )
        
        fname="results/plots/SPOT/CTD_$(varname).png"
    else
        p6, p7 = plot_line_box_violin(wi, sp, su, fa, wi_means, sp_means, su_means, fa_means, xlab, colname, maxdepth)
        l = @layout [ a{0.5w} b
                      c{0.5h} [grid(2, 2)] ]

        fig = plot(p5, p7, p6, p1, p2, p3, p4,
            layout = l,
            fg_legend = :transparent,
            size=(900,900),
            plot_title = "Seasonal CTD $(varname) at SPOT",
        )

        fname="results/plots/SPOT/CTD_$(varname)_$(maxdepth)m.png"
    end

    savefig(fig,fname)

end


function get_ctd_data(maxdepth)

    path = "/home/lee/Dropbox/Development/NPZBD_1D/data/spot_data"

    wi, sp = DataFrame(CSV.File("$(path)/seasonal_ctd/wi_ctd.csv")), DataFrame(CSV.File("$(path)/seasonal_ctd/sp_ctd.csv"))
    su, fa = DataFrame(CSV.File("$(path)/seasonal_ctd/su_ctd.csv")), DataFrame(CSV.File("$(path)/seasonal_ctd/fa_ctd.csv"))
    means = DataFrame(CSV.File("$(path)/seasonal_means/ctd/season_means_ctd.csv"))
    all = DataFrame(CSV.File("$(path)/CTD_Data.csv"))

    if maxdepth != 0
        wi, sp, su, fa, means, all = filter_by_depth([wi, sp, su, fa, means, all], maxdepth)
        return means, wi, sp, su, fa, all
    end

    return means, wi, sp, su, fa, all

end

function get_btl_data()

    path = "/home/lee/Dropbox/Development/NPZBD_1D/data/spot_data"

    wi, sp = DataFrame(CSV.File("$(path)/seasonal/winter.csv")), DataFrame(CSV.File("$(path)/seasonal/spring.csv"))
    su, fa = DataFrame(CSV.File("$(path)/seasonal/summer.csv")), DataFrame(CSV.File("$(path)/seasonal/fall.csv"))
    means = DataFrame(CSV.File("$(path)/seasonal_means/bottle/seasonal_means.csv"))

    return means, wi, sp, su, fa, all

end


function filter_by_depth(data_arr, maxdepth)

    for ds in data_arr
        filter!(row -> row.CTDDEPTH < maxdepth, ds)
    end

    return data_arr

end

# max_depth=60
# var_name="Fluorescence"
# col_name="CTDFLUOR"
# x_lab=L"\mu g/L"

# max_depth=100
# var_name="Attenuation"
# col_name="CTDATTN"
# x_lab=L"m^{-1}"

# max_depth=400
# var_name="DO"
# col_name="CTDOXY"
# x_lab=L"\mu mol/L"

# plot_ctd_data(max_depth, var_name, col_name, x_lab)
# plot_btl_data



# DCM_sp = dropmissing(DCM_sp, :depth)
# plot(ctd_su[!,"CTDTMP"], -y, alpha=0.7, label=false, grid=false, border=:box, title="SPOT CDT Temperature", ylabel="Depth(m)", xlabel="deg C")

# plot!(wi_means[!, "CTDTMP"], -wi_means[!, "CTDDEPTH"],alpha=0.7, label=false, grid=false, border=:box)

# plot!(sp_means[!, "CTDTMP"], -sp_means[!, "CTDDEPTH"],alpha=0.7, label=false, grid=false, border=:box)

# plot!(ctd_fa[!, "CTDTMP"], -ctd_fa[!, "CTDDEPTH"],alpha=0.7, label=false, grid=false, border=:box)




