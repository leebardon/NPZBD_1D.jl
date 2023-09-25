
# using NCDatasets
# using Plots, Colors, LaTeXStrings
# using DataFrames

# include("utils/utils.jl")
# include("utils/save_utils.jl")


function plot_biomasses(fsaven, season_num)

    H = 890
    zc = get_zc(H)

    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    season_num == 1 ? season = "Winter" : season = "Summer"

    ep = get_endpoints(["n", "p", "z", "b", "d", "o"], ds)
    n, p, z, b, d, o = ep[1], ep[2], ep[3], ep[4], ep[5], ep[6]
    
    f1, dir1 = plot_stacked(p, b, z, d, n, o, zc, filename)
    f2, dir2 = plot_individual(b, d, z, p, zc, season, filename)

    savefig(f1,"$(dir1)/$(filename).png")
    savefig(f2,"$(dir2)/$(filename).png")

end


function plot_stacked(p, b, z, d, n, o, zc, filename)

    parent_folder = "results/plots/bmass_stacked/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-500.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    
    p1 = plot(sum(p, dims = 2), -zc, lw=ls, lc="darkgreen", grid=false, xrotation=45, label=" Total P", 
    ylimits=yl, alpha=ab, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mm/m^3")
    plot!(sum(b, dims = 2), -zc, lw=ls, lc="skyblue", label=" Total B", alpha=ab, labelfontsize=lfs, legend=lg) 
    plot!(sum(z, dims = 2), -zc, lw=ls, lc="maroon", label=" Total Z", alpha=ab, labelfontsize=lfs, legend=lg) 

    p2 = plot(sum(d, dims = 2), -zc, lw=ls, lc="orange", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")

    p3 = plot(sum(n, dims = 2), -zc, lw=ls, lc="grey", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")

    p4 = plot(sum(o, dims = 2), -zc, lw=ls, lc="pink", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, yformatter=Returns(""))

    f1 = plot(p1, p2, p3, p4,
        layout = [1 1 1 1],
        size=(700,450),
        fg_legend = :transparent,
        title = ["Biomass" "D" "N" "O"]
    )

    return f1, dir
end


function plot_individual(b, d, z, p, zc, season, filename)

    parent_folder = "results/plots/bmass_individual/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-500.0, 0)
    sizes = get_size([b, d, z, p])
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    nb, nd, nz, np = sizes[1], sizes[2], sizes[3], sizes[4]
    tfs = 9

    # Plot OM
    tls = ["  Organic Matter", "Z Biomass", "B Biomass", "P Biomass"]
    p1 = plot(d[:,1], -zc, lc=dcols[1], lw=ls, grid=false,  label=" D1", xrotation=45, ylimits=yl, title=tls[1], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mm/m^3")
    if nd > 1
        for i in 2:nd
            plot!(d[:,i], -zc, lc=dcols[i], lw=ls, label=" D$i", alpha=ab, labelfontsize=lfs, legend=lg)
        end
    end
  
    # Plot zoo
    p2 = plot(z[:,1], -zc, lc=zcols[1], grid=false, label=" Z1", lw=ls, xrotation=45, ylimits=yl, title=tls[2], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")
    if nz > 1
        for j in 2:nz
            plot!(z[:,j], -zc, lc=zcols[j], lw=ls, label=" Z$j", alpha=ab, labelfontsize=lfs, legend=lg)
        end
    end

    # Plot bac
    p3 = plot(b[:,1], -zc, lc=bcols[1], grid=false, lw=ls, label=" B1", xrotation=45, ylimits=yl, title=tls[3], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")
    if nb > 1
        for k in 2:nb
            plot!(b[:,k], -zc, lc=bcols[k], lw=ls, label=" B$k", alpha=ab, labelfontsize=lfs , legend=lg)
        end
    end

    # Plot phyto
    p4 = plot(p[:,1], -zc, lc=pcols[1], grid=false, lw=ls, label=" P1", xrotation=45, ylimits=yl, title=tls[4], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mm/m^3")
    if np > 1
        for l in 2:np
            plot!(p[:,l], -zc, lc=pcols[l], lw=ls, label=" P$l", alpha=ab, labelfontsize=lfs, legend=lg)
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


# fsaven = "results/outfiles/Wi100y_230923_17:23_8P6Z13B5D.nc"
# plot_biomasses(fsaven, 1)

