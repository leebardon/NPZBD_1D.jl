using Plots, LatexStrings

function plot_mean_RA(ann_means, wi_means, sp_means, su_means, fa_means)

    clades = ["SAR11_clade", "Flavobacteriales", "SAR324_clade", "SAR202_clade"]
    clrs = ["teal", "red4", "purple4", "orchid"]

    lw1=8
    mks=8

    p1 = plot(ann_means[1:3, clades[1]], -ann_means[1:3, "Depth"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    ylabel="Depth (m)", xlabel="Relative Abundance", title="Annual", xrotation=45)
    plot!(ann_means[1:3, clades[2]], -ann_means[1:3, "Depth"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(ann_means[1:3, clades[3]], -ann_means[1:3, "Depth"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(ann_means[1:3, clades[4]], -ann_means[1:3, "Depth"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(ann_means[1:3, clades[1]], -ann_means[1:3, "Depth"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(ann_means[1:3, clades[2]], -ann_means[1:3, "Depth"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(ann_means[1:3, clades[3]], -ann_means[1:3, "Depth"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(ann_means[1:3, clades[4]], -ann_means[1:3, "Depth"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)

    p2 = plot(wi_means[1:3, clades[1]], -wi_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    ylabel="Depth (m)", xlabel="", title="Winter", xrotation=45)
    plot!(wi_means[1:3, clades[2]], -wi_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(wi_means[1:3, clades[3]], -wi_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(wi_means[1:3, clades[4]], -wi_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(wi_means[1:3, clades[1]], -wi_means[1:3, "DepthBin"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(wi_means[1:3, clades[2]], -wi_means[1:3, "DepthBin"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(wi_means[1:3, clades[3]], -wi_means[1:3, "DepthBin"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(wi_means[1:3, clades[4]], -wi_means[1:3, "DepthBin"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)

    p3 = plot(su_means[1:3, clades[1]], -su_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    yformatter=Returns(""), xlabel="", title="Summer", xrotation=45)
    plot!(su_means[1:3, clades[2]], -su_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(su_means[1:3, clades[3]], -su_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(su_means[1:3, clades[4]], -su_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(su_means[1:3, clades[1]], -su_means[1:3, "DepthBin"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(su_means[1:3, clades[2]], -su_means[1:3, "DepthBin"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(su_means[1:3, clades[3]], -su_means[1:3, "DepthBin"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(su_means[1:3, clades[4]], -su_means[1:3, "DepthBin"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)

    p4 = plot(fa_means[1:3, clades[1]], -fa_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    ylabel="Depth (m)", xlabel="", title="Fall", xrotation=45)
    plot!(fa_means[1:3, clades[2]], -fa_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(fa_means[1:3, clades[3]], -fa_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(fa_means[1:3, clades[4]], -fa_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(fa_means[1:3, clades[1]], -fa_means[1:3, "DepthBin"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(fa_means[1:3, clades[2]], -fa_means[1:3, "DepthBin"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(fa_means[1:3, clades[3]], -fa_means[1:3, "DepthBin"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(fa_means[1:3, clades[4]], -fa_means[1:3, "DepthBin"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)
    
    p5 = plot(sp_means[1:3, clades[1]], -sp_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    yformatter=Returns(""), xlabel="", title="Spring", xrotation=45)
    plot!(sp_means[1:3, clades[2]], -sp_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(sp_means[1:3, clades[3]], -sp_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(sp_means[1:3, clades[4]], -sp_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(sp_means[1:3, clades[1]], -sp_means[1:3, "DepthBin"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(sp_means[1:3, clades[2]], -sp_means[1:3, "DepthBin"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(sp_means[1:3, clades[3]], -sp_means[1:3, "DepthBin"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(sp_means[1:3, clades[4]], -sp_means[1:3, "DepthBin"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)

    l = @layout [ a b ; c d ; e ] 

    f = plot(p2, p3, p4, p5, p1,
    layout = l,
    fg_legend = :transparent,
    size=(600,1000),
    plot_title = "Mean Guild Relative Abundance at SPOT\n",
    )

    savefig(f,"RA.pdf")
end


function plot_monthly_mean_RA()

    clades = ["SAR11_clade", "Flavobacteriales", "SAR324_clade", "SAR202_clade"]
    clrs = ["teal", "red4", "purple4", "orchid"]

    lw1=8
    mks=8

    p1 = plot(ann_means[1:3, clades[1]], -ann_means[1:3, "Depth"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    ylabel="Depth (m)", xlabel="Relative Abundance", title="Annual", xrotation=45)
    plot!(ann_means[1:3, clades[2]], -ann_means[1:3, "Depth"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(ann_means[1:3, clades[3]], -ann_means[1:3, "Depth"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(ann_means[1:3, clades[4]], -ann_means[1:3, "Depth"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(ann_means[1:3, clades[1]], -ann_means[1:3, "Depth"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(ann_means[1:3, clades[2]], -ann_means[1:3, "Depth"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(ann_means[1:3, clades[3]], -ann_means[1:3, "Depth"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(ann_means[1:3, clades[4]], -ann_means[1:3, "Depth"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)

    p2 = plot(wi_means[1:3, clades[1]], -wi_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    ylabel="Depth (m)", xlabel="", title="Winter", xrotation=45)
    plot!(wi_means[1:3, clades[2]], -wi_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(wi_means[1:3, clades[3]], -wi_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(wi_means[1:3, clades[4]], -wi_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(wi_means[1:3, clades[1]], -wi_means[1:3, "DepthBin"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(wi_means[1:3, clades[2]], -wi_means[1:3, "DepthBin"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(wi_means[1:3, clades[3]], -wi_means[1:3, "DepthBin"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(wi_means[1:3, clades[4]], -wi_means[1:3, "DepthBin"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)

    p3 = plot(su_means[1:3, clades[1]], -su_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    yformatter=Returns(""), xlabel="", title="Summer", xrotation=45)
    plot!(su_means[1:3, clades[2]], -su_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(su_means[1:3, clades[3]], -su_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(su_means[1:3, clades[4]], -su_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(su_means[1:3, clades[1]], -su_means[1:3, "DepthBin"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(su_means[1:3, clades[2]], -su_means[1:3, "DepthBin"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(su_means[1:3, clades[3]], -su_means[1:3, "DepthBin"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(su_means[1:3, clades[4]], -su_means[1:3, "DepthBin"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)

    p4 = plot(fa_means[1:3, clades[1]], -fa_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    ylabel="Depth (m)", xlabel="", title="Fall", xrotation=45)
    plot!(fa_means[1:3, clades[2]], -fa_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(fa_means[1:3, clades[3]], -fa_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(fa_means[1:3, clades[4]], -fa_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(fa_means[1:3, clades[1]], -fa_means[1:3, "DepthBin"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(fa_means[1:3, clades[2]], -fa_means[1:3, "DepthBin"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(fa_means[1:3, clades[3]], -fa_means[1:3, "DepthBin"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(fa_means[1:3, clades[4]], -fa_means[1:3, "DepthBin"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)
    
    p5 = plot(sp_means[1:3, clades[1]], -sp_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[1], lw=lw1, alpha=0.3, border=:box, ylims=(-160, 0),
    yformatter=Returns(""), xlabel="", title="Spring", xrotation=45)
    plot!(sp_means[1:3, clades[2]], -sp_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[2], lw=lw1, alpha=0.3)
    plot!(sp_means[1:3, clades[3]], -sp_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[3], lw=lw1, alpha=0.3)
    plot!(sp_means[1:3, clades[4]], -sp_means[1:3, "DepthBin"], grid=false, label=false, lc=clrs[4], lw=lw1, alpha=0.3)
    scatter!(sp_means[1:3, clades[1]], -sp_means[1:3, "DepthBin"], grid=false, color=clrs[1], label="SAR11", markersize=mks, alpha=0.7)
    scatter!(sp_means[1:3, clades[2]], -sp_means[1:3, "DepthBin"], grid=false, color=clrs[2], label="Flavo", markersize=mks, alpha=0.7)
    scatter!(sp_means[1:3, clades[3]], -sp_means[1:3, "DepthBin"], grid=false, color=clrs[3], label="SAR324", markersize=mks, alpha=0.7)
    scatter!(sp_means[1:3, clades[4]], -sp_means[1:3, "DepthBin"], grid=false, color=clrs[4], label="SAR202", markersize=mks, alpha=0.7)

    l = @layout [ a b ; c d ; e ] 

    f = plot(p2, p3, p4, p5, p1,
    layout = l,
    fg_legend = :transparent,
    size=(600,1000),
    plot_title = "Mean Guild Relative Abundance at SPOT\n",
    )

    savefig(f,"RA.pdf")
end


function plot_heatmaps_2014(mar, apr)

    #NOTE observations are at discrete depths - for heatmaps, we need to interpolate between depths 

    

end