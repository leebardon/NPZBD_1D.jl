using NCDatasets, Plots, LaTeXStrings

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")

function create_gif(P)
    # gif of first 2 years, plotted once per day
    # gr(fmt = :png)
    d = 200
    zc = get_zc(d)
    zd = Int(d/10)
    Px = P[1:zd,:,1:20:14640]
    ls=8
    tfs=10
    pcols = ["hotpink2", "darkgreen","red4", "cyan4", "gold3", "black", "brown", "wheat2", "mediumpurple3", "darkseagreen" ]

    xl1 = maximum(Px[:,1,:]) ; xl1 < 1 ? xl1 = 1 : xl1 = xl1
    xl2 = maximum(Px[:,2,:]) ; xl2 < 1 ? xl2 = 1 : xl2 = xl2
    xl3 = maximum(Px[:,3,:]) ; xl3 < 1 ? xl3 = 1 : xl3 = xl3
    xl4 = maximum(Px[:,4,:]) ; xl4 < 1 ? xl4 = 1 : xl4 = xl4
    xl5 = maximum(Px[:,5,:]) ; xl5 < 1 ? xl5 = 1 : xl5 = xl5
    xl6 = maximum(Px[:,6,:]) ; xl6 < 1 ? xl6 = 1 : xl6 = xl6

    anim = @animate for i in 1:732
        @views p1, p2, p3, p4, p5, p6 = Px[:,1,i], Px[:,2,i], Px[:,3,i], Px[:,4,i], Px[:,5,i], Px[:,6,i]
        pl1 = plot(p1, -zc, label="", lw=ls, title="P1", titlefontsize=tfs, grid=false, lc=pcols[1], xlims=(0, xl1), xrotation=45, 
        ylabel="Depth (m)", xlabel="")
        pl2 = plot(p2, -zc, label="", lw=ls, title="P2", titlefontsize=tfs, grid=false, lc=pcols[2], xlims=(0, xl2), xrotation=45, 
        yformatter=Returns(""), xlabel="")
        pl3 = plot(p3, -zc, label="", lw=ls, title="P3", titlefontsize=tfs, grid=false, lc=pcols[3], xlims=(0, xl3), xrotation=45, 
        yformatter=Returns(""), xlabel="")
        pl4 = plot(p4, -zc, label="", lw=ls, title="P4", titlefontsize=tfs, grid=false, lc=pcols[4], xlims=(0, xl4), xrotation=45, 
        ylabel="Depth (m)", xlabel=L" mmol ~N/m^3")
        pl5 = plot(p5, -zc, label="", lw=ls, title="P5", titlefontsize=tfs, grid=false, lc=pcols[5], xlims=(0, xl5), xrotation=45, 
        yformatter=Returns(""), xlabel=L" mmol ~N/m^3")
        pl6 = plot(p6, -zc, label="", lw=ls, title="P6", titlefontsize=tfs, grid=false, lc=pcols[6], xlims=(0, xl6), xrotation=45, 
        legend = :bottomright, labelfontsize=10, yformatter=Returns(""), xlabel=L" mmol ~N/m^3")
        plot!(p6, -zc, label=" Day $i", lw=ls, lc=pcols[6], alpha=0, xlims=(0, xl6), legend = :bottomright, legendfontsize=8, 
        yformatter=Returns(""), xlabel=L" mmol ~N/m^3")

        layout = (2,1)
        plot(pl1, pl2, pl3, pl4, pl5, pl6, size=(400, 450), fg_legend = :transparent)
    end

    dir = "results/plots/animations"
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")

    gif(anim, "$(dir)/P_$(filename).gif", fps=25)

end



# fsaven="results/outfiles/Wi50y_231213_14:58_6P3Z13B8D.nc"
# fsaven="results/outfiles/Su50y_231213_15:32_6P3Z13B8D.nc"

# fsaven="results/outfiles/Wi50y_231214_20:10_6P3Z13B8D.nc"

# ds = NCDataset(fsaven)
# P = ds["p"][:,:,:]

# create_gif(P)
