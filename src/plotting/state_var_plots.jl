
using NCDatasets
using Plots, Colors, LaTeXStrings
using DataFrames, Statistics

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")


function plot_state_vars(fsaven)

    H = 890
    zc = get_zc(H)

    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    ds.attrib["Season"] == "winter" ? season_num == 1 : season_num == 2

    if ds["pulse"][:] == 1
        N, P, Z, B, D, O = get_endpoints(["n", "p", "z", "b", "d", "o"], ds)
    else
        N, P, Z, B, D, O = mean_over_time(["n", "p", "z", "b", "d", "o"], ds, season_num)
    end
    
    f1, dir1 = plot_stacked(N, P, Z, B, D, O, zc, filename)
    savefig(f1,"$(dir1)/$(filename).png")

    f2, dir2 = plot_individual(P, Z, B, D, zc, season, filename)
    savefig(f2,"$(dir2)/$(filename).png")

end


function plot_stacked(N, P, Z, B, D, O, zc, filename)

    parent_folder = "results/plots/bmass_stacked/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-890.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls=7
    ab=1.0
    lfs = 8
    
    p1 = plot(sum(P, dims=2), -zc, lw=ls, lc="darkgreen", grid=false, xrotation=45, label=" Total P", 
    ylimits=yl, alpha=ab, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mmol N/m^3")
    plot!(sum(B, dims=2), -zc, lw=ls, lc="skyblue", label=" Total B", alpha=ab, labelfontsize=lfs, legend=lg) 
    plot!(sum(Z, dims=2), -zc, lw=ls, lc="maroon", label=" Total Z", alpha=ab, labelfontsize=lfs, legend=lg) 

    p2 = plot(sum(D, dims=2), -zc, lw=ls, lc="orange", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")

    p3 = plot(sum(N, dims=2), -zc, lw=ls, lc="grey", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")

    p4 = plot(O[:,1], -zc, lw=ls, lc="pink", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol N/m^3")

    f1 = plot(p1, p2, p3, p4,
        layout = [1 1 1 1],
        size=(600,350),
        fg_legend = :transparent,
        title = ["Biomass" "D" "N" "O"]
    )

    return f1, dir
end


function plot_individual(P, Z, B, D, zc, season, filename)

    parent_folder = "results/plots/bmass_individual/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-890.0, 0)
    sizes = get_size([B, D, Z, P])
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    nb, nd, nz, np = sizes[1], sizes[2], sizes[3], sizes[4]
    tfs = 9
    ls=7
    ab=1.0
    lfs = 8

    # Plot OM
    tls = ["  Organic Matter", "Z Biomass", "B Biomass", "P Biomass"]
    p1 = plot(D[:,1], -zc, lc=dcols[1], lw=ls, grid=false,  label=" D1", xrotation=45, ylimits=yl, title=tls[1], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mmol N/m^3")
    if nd > 1
        for i in 2:nd
            plot!(D[:,i], -zc, lc=dcols[i], lw=ls, label=" D$i", alpha=ab, labelfontsize=lfs, legend=lg)
        end
    end
  
    # Plot zoo
    p2 = plot(Z[:,1], -zc, lc=zcols[1], grid=false, label=" Z1", lw=ls, xrotation=45, ylimits=yl, title=tls[2], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")
    if nz > 1
        if nz > 6
            for j in 2:nz
                plot!(Z[:,j], -zc, lw=ls, label=" Z$j", alpha=ab, labelfontsize=lfs, legend=lg)
            end
        else
            for j in 2:nz
                plot!(Z[:,j], -zc, lc=zcols[j], lw=ls, label=" Z$j", alpha=ab, labelfontsize=lfs, legend=lg)
            end
        end
    end

    # Plot bac
    bcols = ["red3", "darkorange", "indigo", "lightgreen", "olivedrab4", "darkgreen", "cyan2", "dodgerblue",
    "gray85", "gray67", "gray48",  "gray30", "black"] #uncomment for darwin runs
    p3 = plot(B[:,1], -zc, lc=bcols[1], grid=false, lw=ls, linestyle=:dot, label=" B1", xrotation=45, ylimits=yl, title=tls[3], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")
    if nb > 1
        for k in 2:nb
            if k == 2 || k == 3
                plot!(B[:,k], -zc, lc=bcols[k], lw=ls, linestyle=:dot, label=" B$k", alpha=ab, labelfontsize=lfs , legend=lg)
            elseif k == 4 
                plot!(B[:,k], -zc, lc=bcols[k], lw=ls, label=" B$(k)", alpha=ab, labelfontsize=lfs, legend=lg)
            elseif 5 <= k <= 8
                plot!(B[:,k], -zc, lc=bcols[k], lw=ls, label=" B$k", alpha=ab, labelfontsize=lfs, legend=lg)
            elseif k == 9
                plot!(B[:,k], -zc, lc=bcols[k], lw=ls, label=" B$(k)", alpha=ab, labelfontsize=lfs, legend=lg)
            else
                plot!(B[:,k], -zc, lc=bcols[k], lw=ls, label=" B$k", alpha=ab, labelfontsize=lfs, legend=lg)
            end
        end
    end

    # Plot phyto
    p4 = plot(P[:,1], -zc, lc=pcols[1], grid=false, lw=ls, label=" P1", xrotation=45, ylimits=yl, title=tls[4], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")
    if np > 1
        for l in 2:np
            plot!(P[:,l], -zc, lc=pcols[l], lw=ls, label=" P$l", alpha=ab, labelfontsize=lfs, legend=lg)
        end
    end

    f2 = plot(p1, p2, p3, p4,
            layout = [1 1 1 1],
            fg_legend = :transparent,
            size=(600,350),
            # plot_title = season,
        )
    
    return f2, dir

end

##------------------------------------------------------------------------------------------------------------


# fsaven = "results/outfiles/Wi100y_231011_20:23_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi100y_231011_23:28_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi30y_231012_22:36_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi2y_231013_14:09_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi100y_231017_01:23_10P3Z21B9D.nc"
# fsaven = "results/outfiles/Wi100y_231017_12:01_10P3Z21B9D.nc"

# fsaven = "results/outfiles/Wi50y_231118_10:40_6P3Z13B8D.nc" #meso pulse 50yrs
# fsaven = "results/outfiles/Wi50y_231118_14:34_6P3Z13B8D.nc" #meso steady 50yrs
# fsaven = "results/outfiles/Su50y_231118_14:08_6P3Z13B8D.nc" #oli pulse 50yrs
# fsaven = "results/outfiles/Su50y_231118_15:02_6P3Z13B8D.nc" #oli steady 50yrs

# fsaven = "results/outfiles/Wi100y_231119_13:54_6P3Z13B8D.nc" # meso pulse 100yrs
# fsaven = "results/outfiles/Wi100y_231119_15:34_6P3Z13B8D.nc" # meso steady 100yrs
# fsaven = "results/outfiles/Su100y_231119_14:25_6P3Z13B8D.nc" #oligo pulse 100yrs
# fsaven = "results/outfiles/Su100y_231119_16:05_6P3Z13B8D.nc" #oligo steady 100yrs

# plot_state_vars_dar(fsaven, 1, "")
# plot_state_vars_dar(fsaven, 1, "(steady state)")
# plot_state_vars_dar(fsaven, 2, "")
# plot_state_vars_dar(fsaven, 2, "(steady state)")


# fsaven = "results/outfiles/Wi50y_231120_09:29_6P2Z10B10D.nc" # meso pulse 100yrs - lability test
# plot_state_vars(fsaven, 1)
# plot_state_vars_lab(fsaven, 1)