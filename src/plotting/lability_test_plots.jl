
function plot_state_vars_lab(fsaven, season_num)

    H = 890
    zc = get_zc(H)

    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    season_num == 1 ? season = "Mesotrophic" : season = "Oligotrophic"

    if ds["pulse"][:] == 1
        N, P, Z, B, D, O = get_endpoints(["n", "p", "z", "b", "d", "o"], ds)
    else
        N, P, Z, B, D, O = mean_over_time(["n", "p", "z", "b", "d", "o"], ds, season)
    end
    
    f1, dir1 = plot_stacked_lab(N, P, Z, B, D, O, zc, filename)
    println("saving to: $(dir1)/$(filename)_lab.png")
    savefig(f1,"$(dir1)/$(filename)_lab.png")

    f2, dir2 = plot_individual_lab(P, Z, B, D, N, zc, season, filename)
    println("saving to: $(dir2)/$(filename)_lab.png")
    savefig(f2,"$(dir2)/$(filename)_lab.png")

end


function plot_stacked_lab(N, P, Z, B, D, O, zc, filename)

    parent_folder = "results/plots/bmass_stacked/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-890.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls=7
    ab=1.0
    lfs = 8
    
    p1 = plot(sum(P, dims=2), -zc, lw=ls, lc="darkgreen", grid=false, xrotation=45, label=" Total P", 
    ylimits=yl, alpha=ab, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L"mmol N/m^3")
    plot!(sum(B, dims=2), -zc, lw=ls, lc="skyblue", label=" Total B", alpha=ab, labelfontsize=lfs, legend=lg) 
    plot!(sum(Z, dims=2), -zc, lw=ls, lc="maroon", label=" Total Z", alpha=ab, labelfontsize=lfs, legend=lg) 

    p2 = plot(sum(D, dims=2), -zc, lw=ls, lc="orange", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L"mmol N/m^3")

    p3 = plot(sum(N, dims=2), -zc, lw=ls, lc="grey", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L"mmol N/m^3")

    p4 = plot(O[:,1], -zc, lw=ls, lc="pink", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, yformatter=Returns(""), xlabel=L"mmol N/m^3")

    f1 = plot(p1, p2, p3, p4,
        layout = [1 1 1 1],
        size=(600,350),
        fg_legend = :transparent,
        title = ["Biomass" "D" "N" "O"]
    )

    return f1, dir
end



function plot_individual_lab(P, Z, B, D, N, zc, season, filename)

    parent_folder = "results/plots/bmass_individual/"
    dir = check_subfolder_exists(filename, parent_folder)

    P, Z, B = set_extinct_to_zero(P), set_extinct_to_zero(Z), set_extinct_to_zero(B)

    yl=(-890.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls=7
    ab=0.6
    lfs = 6
    ls2=5
    
    dcols = ["teal", "lightblue1", "azure4", "red4", "black", "magenta2", "thistle", "seagreen4", "goldenrod4", "grey"]
    bcols = ["teal", "lightblue1", "azure4", "red4", "black", "magenta2", "thistle", "seagreen4", "goldenrod4", "grey"]

    p1 = plot(D[:,1], -zc, lc=dcols[1], lw=ls2, linestyle=:dot, grid=false,  label=" Most Lab", xrotation=45, ylimits=yl, title="L to SL POM", 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mmol N/m^3")
    plot!(D[:,2], -zc, lc=dcols[2], lw=ls2, linestyle=:dot, label=" ", labelfontsize=lfs, legend=lg)
    plot!(D[:,3], -zc, lc=dcols[3], lw=ls2, linestyle=:dot, label=" ", labelfontsize=lfs, legend=lg)
    plot!(D[:,4], -zc, lc=dcols[4], lw=ls2, linestyle=:dot, label=" ", labelfontsize=lfs, legend=lg)
    plot!(D[:,5], -zc, lc=dcols[5], lw=ls2, linestyle=:dot, label=" Least Lab", labelfontsize=lfs, legend=lg)

    p2 = plot(D[:,6], -zc, lc=dcols[6], lw=ls2, linestyle=:dot, grid=false,  label=" Least Rec", xrotation=45, ylimits=yl, title="Semi-Rec. to Rec. POM", 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mmol N/m^3")
    plot!(D[:,7], -zc, lc=dcols[7], lw=ls2, linestyle=:dot, label=" ", labelfontsize=lfs, legend=lg)
    plot!(D[:,8], -zc, lc=dcols[8], lw=ls2, linestyle=:dot, label=" ", labelfontsize=lfs, legend=lg)
    plot!(D[:,9], -zc, lc=dcols[9], lw=ls2, linestyle=:dot, label=" ", labelfontsize=lfs, legend=lg)
    plot!(D[:,10], -zc, lc=dcols[10], lw=ls2, linestyle=:dot, label=" Most Rec", labelfontsize=lfs, legend=lg)

    p3 = plot(B[:,1], -zc, lc=bcols[1], lw=ls, grid=false,  label=" B1", xrotation=45, ylimits=yl, title="Lab to Semi-Lab POM Cons", 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")
    plot!(B[:,2], -zc, lc=bcols[2], lw=ls, label=" B2", labelfontsize=lfs , legend=lg)
    plot!(B[:,3], -zc, lc=bcols[3], lw=ls, label=" B3", labelfontsize=lfs , legend=lg)
    plot!(B[:,4], -zc, lc=bcols[4], lw=ls, label=" B4", labelfontsize=lfs , legend=lg)
    plot!(B[:,5], -zc, lc=bcols[5], lw=ls, label=" B5", labelfontsize=lfs , legend=lg)

    p4 = plot(B[:,6], -zc, lc=bcols[6], lw=ls, grid=false,  label=" B6", xrotation=45, ylimits=yl, title="Semi-Rec. to Rec. POM Cons", 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")
    plot!(B[:,7], -zc, lc=bcols[7], lw=ls, label=" B7", labelfontsize=lfs , legend=lg)
    plot!(B[:,8], -zc, lc=bcols[8], lw=ls, label=" B8", labelfontsize=lfs , legend=lg)
    plot!(B[:,9], -zc, lc=bcols[9], lw=ls, label=" B9", labelfontsize=lfs , legend=lg)
    plot!(B[:,10], -zc, lc=bcols[10], lw=ls, label=" B10", labelfontsize=lfs , legend=lg)


    f2 = plot(p1, p2, p3, p4,
            layout = (1, 4),
            fg_legend = :transparent,
            size=(600,350),
            # plot_title = "$season $type",
        )
    
    return f2, dir

end