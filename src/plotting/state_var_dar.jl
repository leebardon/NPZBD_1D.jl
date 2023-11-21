function plot_state_vars_dar(fsaven, season_num, type)

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
    
    f1, dir1 = plot_stacked_dar(N, P, Z, B, D, O, zc, filename)
    println("saving to: $(dir1)/$(filename)_dar.png")
    savefig(f1,"$(dir1)/$(filename)_dar.png")

    f2, dir2 = plot_individual_dar(P, Z, B, D, N, zc, season, filename, type)
    println("saving to: $(dir2)/$(filename)_dar.png")
    savefig(f2,"$(dir2)/$(filename)_dar.png")

end


function plot_stacked_dar(N, P, Z, B, D, O, zc, filename)

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



function plot_individual_dar(P, Z, B, D, N, zc, season, filename, type)

    parent_folder = "results/plots/bmass_individual/"
    dir = check_subfolder_exists(filename, parent_folder)

    P, Z, B = set_extinct_to_zero(P), set_extinct_to_zero(Z), set_extinct_to_zero(B)

    yl=(-890.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls=7
    ab=0.6
    lfs = 8
    ls2=5
    
    dcols = ["teal", "lightblue1", "azure4", "red4", "black"]
    bcols = ["teal", "lightblue1", "azure4", "red4", "pink", "black","grey"]

    tls = ["POM", "POM Consumers", "DOM", "DOM Consumers"]
    lab_om=sum([D[:,4], D[:,5], D[:,6]])
    lab_cop=sum([B[:,4], B[:,5], B[:,6]])
    lab_oli=sum([B[:,9], B[:,10], B[:,11]])

    p1 = plot(D[:,1], -zc, lc=dcols[1], lw=ls2, linestyle=:dot, grid=false,  label=" Labile", xrotation=45, ylimits=yl, title=tls[1], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mmol N/m^3")
    plot!(D[:,2], -zc, lc=dcols[2], lw=ls2, linestyle=:dot, label=" ", labelfontsize=lfs, legend=lg)
    plot!(D[:,3], -zc, lc=dcols[3], lw=ls2, linestyle=:dot, label=" Recalc.", labelfontsize=lfs, legend=lg)

    p2 = plot(B[:,1], -zc, lc=bcols[1], lw=ls, grid=false,  label=" Copio.", xrotation=45, ylimits=yl, title=tls[2], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")
    plot!(B[:,2], -zc, lc=bcols[2], lw=ls, label=" ", labelfontsize=lfs , legend=lg)
    plot!(B[:,3], -zc, lc=bcols[3], lw=ls, label=" Oligo.", labelfontsize=lfs , legend=lg)

    p3 = plot(lab_om, -zc, lc=dcols[4], lw=ls2, linestyle=:dot, grid=false, label=" Labile", xrotation=45, ylimits=yl, title=tls[3], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")
    plot!(D[:,7], -zc, lc=dcols[5], lw=ls2, linestyle=:dot,label=" Recalc.", labelfontsize=lfs, legend=lg)

    p4 = plot(lab_cop, -zc, lc=bcols[4], grid=false, lw=ls, label=" Copio.", xrotation=45, ylimits=yl, title=tls[4], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3", alpha=ab)
    plot!(lab_oli, -zc, lc=bcols[5], lw=ls, label=" Oligo.",labelfontsize=lfs , legend=lg, alpha=ab)
    plot!(B[:,7], -zc, lc=bcols[6], lw=ls, label=" Copio.",labelfontsize=lfs , legend=lg, alpha=ab)
    plot!(B[:,12], -zc, lc=bcols[7], lw=ls, label=" Oligo.", labelfontsize=lfs , legend=lg, alpha=ab)


    f2 = plot(p1, p2, p3, p4,
            layout = (1, 4),
            fg_legend = :transparent,
            size=(600,350),
            # plot_title = "$season $type",
        )
    
    return f2, dir

end