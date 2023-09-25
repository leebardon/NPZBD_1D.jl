using NCDatasets
using Plots, Colors, LaTeXStrings
using DataFrames

bc= ["cyan3", "darkorange", "indigo", "coral4", "lightcyan4", "magenta2", "orange4", "seagreen4",
"darkkhaki", "purple", "crimson",  "azure4", "turquoise1"]
dc= ["blue3", "black", "maroon", "navy", "brown4"]
pc = ["olivedrab3", "darkgreen","red4", "cyan4", "purple", "black", "hotpink2", "wheat2" ]
nc = ["blue2"]
ab=0.8
ab_ext=0.8
ls=5
lfs=9


function depth_plots(fsaven, season_num, years, run_type)

    ds = NCDataset(fsaven)
    file = replace(fsaven, "out_$(years)y_" => "", ".nc" => "", "results/outfiles/" => "")

    season_num == 1 ? season = "Winter" : season = "Summer"

    n, p, z, b, d, o = get_endpoints(ds, ["n", "p", "z", "b", "d", "o"])

    H = ds["H"][:]
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]
    
    f1 = plot_stacked(p, b, z, d, n, o, zc, season)
    f2 = plot_individual(n, b, d, z, p, zc, run_type, season)

    combined = plot(f1, f2, 
        layout = (2,1),
        size=(700,1000),
        plot_title = season,
    )

    savefig(f1,"results/plots/stacked/$(file)_$(season)_$(years)y.png")
    savefig(f2,"results/plots/individual/$(file)_$(season)_$(years)y.png")
    savefig(combined,"results/plots/combined/$(file)_$(season)_$(years)y.png")

end

function plot_stacked(p, b, z, d, n, o, zc, season)

    yl=(-400.0, 0)
    p1 = plot(sum(p, dims = 2), -zc, lc="green", grid=false, xrotation=45, label="Total P", ylimits=yl)
    plot!(sum(b, dims = 2), -zc, lc="blue", label="Total B") 
    plot!(sum(z, dims = 2), -zc, lc="black", label="Total Z") 
    p2 = plot(sum(d, dims = 2), -zc, lc="orange", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl)
    p3 = plot(sum(n, dims = 2), -zc, lc="grey", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl)
    p4 = plot(sum(o, dims = 2), -zc, lc="pink", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl)
    f1 = plot(p1, p2, p3, p4,
        linewidth = 2,
        layout = [1 1 1 1],
        size=(700,500),
        fg_legend = :transparent,
        title = ["Biomass" "D" "N" "O2"]
    )

    return f1
end


function plot_individual(n, b, d, z, p, zc, run_type, season)

    yl=(-500.0, 0)
    sizes = get_size([b, d, z, p])
    nb, nd, nz, np = sizes[1], sizes[2], sizes[3], sizes[4]
    tfs = 9

    # Plot OM
    tls = ["  Nutrients", "Z Biomass", "B Biomass", "P Biomass"]
    if run_type == 3
        p1 = plot(d[:,1], -zc, linecolor="yellow2", lw=3, alpha=0.5, label=L" POM", xrotation=45, ylimits=yl, title=tls[1], 
        titlefontsize=tfs, legend_position=:bottomright)
        plot!(d[:,2], -zc, grid=false, linecolor="goldenrod", lw=3, alpha=0.5, label=L" DOM")
        plot!(n[:,1], -zc, grid=false, linecolor="grey20", lw=3, alpha=0.5, label=L" N")
    else
        p1 = plot(d[:,1], -zc, linecolor="rosybrown2", lw=3, alpha=0.5, label=L" D1", xrotation=45, ylimits=yl, title=tls[1], 
        titlefontsize=tfs, legend_position=:bottomright)
        c1, c2 = colorant"tan2", colorant"brown4"
        cmp1 = range(c1, stop=c2, length=nd-1)
        for i in 2:nd
            plot!(d[:,i], -zc, palette=cmp1, grid=false, lw=3, alpha=0.5, label=L" D$i")
        end
    end

    # Plot zoo
    if run_type == 3
        p2 = plot(z[:,1], -zc, linecolor="dodgerblue2", grid=false, label=L" P_{gzr}", lw=3, alpha=0.5, xrotation=45, ylimits=yl, 
        title=tls[2], titlefontsize=tfs, legend_position=:bottomright)
        # plot!(z[:,2], -zc, linecolor="seagreen", lw=3,alpha=0.5,  label=L" B_{gz}")
        plot!(z[:,2], -zc, linecolor="seagreen", lw=3,alpha=0.5,  label=L" B1_{gzr}")
        plot!(z[:,3], -zc, linecolor="red4", lw=3,alpha=0.5,  label=L" B2_{gzr}")
        # plot!(z[:,3], -zc, linecolor="red4", lw=3, alpha=0.5, label=L" BAC_g")
    else
        p2 = plot(z[:,1], -zc, linecolor="salmon1", grid=false, label=L" Z1", lw=3, alpha=0.5,  xrotation=45, ylimits=yl, 
        title=tls[2], titlefontsize=tfs, legend_position=:bottomright)
        c1, c2 = colorant"salmon2", colorant"red4"
        cmp2 = range(c1, stop=c2, length=nz-1)
        for j in 2:nz
            plot!(z[:,j], -zc, palette=cmp2, lw=3, alpha=0.5, label=L" Z$i")
        end
    end

    # Plot bac
    if run_type == 3
        p3 = plot(b[:,1], -zc, linecolor="darkorange", grid=false, lw=3, alpha=0.5, label=L" B1_{POM}", xrotation=45, ylimits=yl, 
        title=tls[3], titlefontsize=tfs, legend_position=:bottomright)
        plot!(b[:,2], -zc, linecolor="lightcyan4", lw=3, alpha=0.5, label=L" B2_{DOM}")
    else
        p3 = plot(b[:,1], -zc, linecolor="plum", grid=false, lw=3,alpha=0.5,  ls=:dot, label="", xrotation=45, ylimits=yl, title=tls[3], titlefontsize=tfs)
        c1, c2 = colorant"plum1", colorant"purple4"
        cmp3 = range(c1, stop=c2, length=nb-1)
        for k in 2:nb
            plot!(b[:,k], -zc, palette=cmp3, lw=3, alpha=0.5, label="")
        end
    end

    # Plot phyto
    p4 = plot(p[:,1], -zc, linecolor="powderblue", grid=false, lw=3, alpha=0.5, label=L" P1", xrotation=45, ylimits=yl, 
    title=tls[4], titlefontsize=tfs, legend_position=:bottomright)
    # plot!(p[:,2], -zc, linecolor="maroon", lw=3, alpha=0.5, label=L" P2")
    # plot!(p[:,2], -zc, linecolor="hotpink2", lw=3, alpha=0.5, label=L" P3")
    # plot!(p[:,2], -zc, linecolor="rosybrown2", lw=3, alpha=0.5, label=L" P4")
    
    # c1, c2 = colorant"pink", colorant"maroon"
    # cmp4 = range(c1, stop=c2, length=np-2)
    # for l in 2:np
    #     if l == np
    #         plot!(p[:,l], -zc, grid=false, linecolor="black", lw=3, label=" Rate opt")
    #     else
    #         plot!(p[:,l], -zc, palette=cmp4, grid=false, lw=3, label="")
    #     end
    # end

    f2 = plot(p1, p2, p3, p4,
            layout = [1 1 1 1],
            fg_legend = :transparent,
            size=(700,450),
            plot_title = season,
        )
    
    return f2

end

# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------



# default(show = true)
# outfile = "results/outfiles/out_100y_20230906_1413.nc"
# depth_plots(outfile, 2, 100, 3)
