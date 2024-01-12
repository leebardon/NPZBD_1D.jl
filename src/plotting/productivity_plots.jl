# using Plots, LatexStrings

function plot_seasonal_BP(BP_model, spot_leu, spot_thy, ds, fsaven, season_num)

    H = 890
    zc = get_zc(H)

    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    parent_folder = "results/plots/productivity/"
    dir = check_subfolder_exists(filename, parent_folder)

    season_num == 1 ? season = "Meso." : season = "Oligo."
    leu, thy = dropmissing(spot_leu, disallowmissing=true), dropmissing(spot_thy, disallowmissing=true)

    tfs = 9
    ls = 5
    ls2 = 7
    ab = 0.5
    lfs = 7
    ls3 = 3
    xtfs = 8
    lg = :bottomright

    bcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
    tls = ["\nPOM Consumers (model)", "\nDOM Consumers (model)", "\nTotal", "SPOT (mean_Leu)", "SPOT (mean_Thy)"]

    p1 = plot(BP_model[:, 1], -zc, lw=ls2, lc=bcols[1], label=" B1", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L"cells/mL/day", 
        xrotation=45, title =tls[1], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab)
        plot!(BP_model[:, 2], -zc, lw=ls2, lc=bcols[2], label=" B2", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 3], -zc, lw=ls2, lc=bcols[3], label=" B3", legendfontsize=lfs, alpha=ab)

    p2 = plot(BP_model[:, 4], -zc, lw=ls2, lc=bcols[4], label=" B4c", legendfontsize=lfs, yformatter=Returns(""), xlabel=L"cells/mL/day" ,
        xrotation=45, title =tls[2], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 9], -zc, lw=ls2, lc=bcols[9], label=" B9o", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 5], -zc, lw=ls2, lc=bcols[5], label=" B5c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 10], -zc, lw=ls2, lc=bcols[10], label=" B10o", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 6], -zc, lw=ls2, lc=bcols[6], label=" B6c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 11], -zc, lw=ls2, lc=bcols[11], label=" B11o", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 7], -zc, lw=ls2, lc=bcols[7], label=" B7c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 12], -zc, lw=ls2, lc=bcols[12], label=" B12o", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 8], -zc, lw=ls2, lc=bcols[8], label=" B8c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 13], -zc, lw=ls2, lc=bcols[13], label=" B13o", legendfontsize=lfs, alpha=ab)

    p3 = plot(sum(BP_model[:, 1:3], dims=2), -zc, lw=ls2, lc="cyan", label=" POM", legendfontsize=lfs, yformatter=Returns(""),  xlabel=L"cells/mL/day" ,
        xrotation=45, title =tls[3], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab)
        plot!(sum(BP_model[:, 4:13], dims=2), -zc, lw=ls2, lc="purple", label=" DOM", legendfontsize=lfs, alpha=ab)
        plot!(sum(BP_model[:, :], dims=2), -zc, lw=ls, lc="red", label=" Total", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        scatter!(leu.mean_Leu, -leu.depth,  markersize=7, markercolor="yellow2", markerstrokecolor=:black, markershape=:diamond, 
        label=" SPOT")
        # plot!(p4, ylims=yl, yformatter=Returns(""))


    #------------------------------------------------------------------------------------------------------------------------------
    # leu, thy = dropmissing(spot_leu, disallowmissing=true), dropmissing(spot_thy, disallowmissing=true)
    # yl=(-890.0, 5)

    # p4 = plot(1, type="n", xlab="", ylab="", ylim=c(-890.0, 5), grid=false, xlabel=L"cells/mL/day", yformatter=:none, 
    #     xrotation=45, title =tls[4], titlefontsize=tfs, border=:box, legend=false, xtickfontsize=xtfs)
    #     plot!(leu.mean_Leu, -leu.depth, xerr=leu.sd_Leu; marker=(:circle,5), lw=ls3, lc="red", linestyle=:dot, alpha=ab)
        # plot!(thy.mean_Thy, -thy.depth, xerr=thy.sd_Thy; marker=(:circle,6), lw=ls, lc="grey", alpha=ab, linestyle=:dot)

    # p4 = plot(leu.mean_Leu, -leu.depth, xerr=leu.sd_Leu; marker=(:circle,6), lw=2, lc="red", grid=false, xlabel=L"cells/mL/day",
    #     xrotation=45, title =tls[4], titlefontsize=tfs, border=:box, widen=true,legend=false, xtickfontsize=xtfs, linestyle=:dot, alpha=0.3)
        # plot!(p4, ylims=yl, yformatter=Returns(""))

    # p5 = plot(thy.mean_Thy, -thy.depth, xerr=thy.sd_Thy; marker=(:circle,6), lw=ls, lc="red", grid=false, yformatter=Returns(""), xlabel=L"cells/mL/day",
    #     xrotation=45, title =tls[4], titlefontsize=tfs, border=:box, legend=false, xtickfontsize=xtfs, alpha=ab, linestyle=:dot)

    f = plot(p1, p2, p3,
    layout = [1 1 1],
    fg_legend = :transparent,
    size=(600,550),
    plot_title = "Bacterial Productivity ($season)\n",
    )

    return f, dir, filename 

end