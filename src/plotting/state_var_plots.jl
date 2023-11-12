
using NCDatasets
using Plots, Colors, LaTeXStrings
using DataFrames, Statistics

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")


function plot_state_vars(fsaven, season_num)

    H = 890
    zc = get_zc(H)

    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    season_num == 1 ? season = "Winter" : season = "Summer"

    if ds["pulse"][:] == 1
        N, P, Z, B, D, O = get_endpoints(["n", "p", "z", "b", "d", "o"], ds)
    else
        N, P, Z, B, D, O = mean_over_time(["n", "p", "z", "b", "d", "o"], ds, season)
    end
    
    f1, dir1 = plot_stacked(N, P, Z, B, D, O, zc, filename)
    f2, dir2 = plot_individual(P, Z, B, D, zc, season, filename)

    savefig(f1,"$(dir1)/$(filename).png")
    savefig(f2,"$(dir2)/$(filename).png")

end


function plot_stacked(N, P, Z, B, D, O, zc, filename)

    parent_folder = "results/plots/bmass_stacked/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-500.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    
    p1 = plot(sum(P, dims=2), -zc, lw=ls, lc="darkgreen", grid=false, xrotation=45, label=" Total P", 
    ylimits=yl, alpha=ab, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mm/m^3")
    plot!(sum(B, dims=2), -zc, lw=ls, lc="skyblue", label=" Total B", alpha=ab, labelfontsize=lfs, legend=lg) 
    plot!(sum(Z, dims=2), -zc, lw=ls, lc="maroon", label=" Total Z", alpha=ab, labelfontsize=lfs, legend=lg) 

    p2 = plot(sum(D, dims=2), -zc, lw=ls, lc="orange", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")

    p3 = plot(sum(N, dims=2), -zc, lw=ls, lc="grey", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")

    p4 = plot(O[:,1], -zc, lw=ls, lc="pink", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, yformatter=Returns(""), xlabel=L" mm/m^3")

    f1 = plot(p1, p2, p3, p4,
        layout = [1 1 1 1],
        size=(700,450),
        fg_legend = :transparent,
        title = ["Biomass" "D" "N" "O"]
    )

    return f1, dir
end


function plot_individual(P, Z, B, D, zc, season, filename)

    parent_folder = "results/plots/bmass_individual/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-500.0, 0)
    sizes = get_size([B, D, Z, P])
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    nb, nd, nz, np = sizes[1], sizes[2], sizes[3], sizes[4]
    tfs = 9

    # Plot OM
    tls = ["  Organic Matter", "Z Biomass", "B Biomass", "P Biomass"]
    p1 = plot(D[:,1], -zc, lc=dcols[1], lw=ls, grid=false,  label=" D1", xrotation=45, ylimits=yl, title=tls[1], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mm/m^3")
    if nd > 1
        for i in 2:nd
            plot!(D[:,i], -zc, lc=dcols[i], lw=ls, label=" D$i", alpha=ab, labelfontsize=lfs, legend=lg)
        end
    end
  
    # Plot zoo
    p2 = plot(Z[:,1], -zc, lc=zcols[1], grid=false, label=" Z1", lw=ls, xrotation=45, ylimits=yl, title=tls[2], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")
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
    p3 = plot(B[:,1], -zc, lc=bcols[1], grid=false, lw=ls, label=" B1", xrotation=45, ylimits=yl, title=tls[3], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")
    if nb > 1
        for k in 2:nb
            plot!(B[:,k], -zc, lc=bcols[k], lw=ls, label=" B$k", alpha=ab, labelfontsize=lfs , legend=lg)
        end
    end

    # Plot phyto
    p4 = plot(P[:,1], -zc, lc=pcols[1], grid=false, lw=ls, label=" P1", xrotation=45, ylimits=yl, title=tls[4], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")
    if np > 1
        for l in 2:np
            plot!(P[:,l], -zc, lc=pcols[l], lw=ls, label=" P$l", alpha=ab, labelfontsize=lfs, legend=lg)
        end
    end

    f2 = plot(p1, p2, p3, p4,
            layout = [1 1 1 1],
            fg_legend = :transparent,
            size=(700,450),
            plot_title = season,
        )
    
    return f2, dir

end


# fsaven = "results/outfiles/Wi100y_231011_20:23_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi100y_231011_23:28_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi30y_231012_22:36_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi2y_231013_14:09_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi100y_231017_01:23_10P3Z21B9D.nc"
# fsaven = "results/outfiles/Wi100y_231017_12:01_10P3Z21B9D.nc"
fsaven = "results/outfiles/Wi2y_231112_13:32_8P6Z13B5D.nc"
plot_state_vars(fsaven, 1)

