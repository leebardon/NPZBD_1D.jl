
using NCDatasets
using Plots, Colors, LaTeXStrings
using DataFrames, Statistics

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/rstar.jl")


function plot_state_vars_dar(fsaven, season_num)

    H = 890
    zc = get_zc(H)

    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    ds.attrib["Season"] == "winter" ? season_num == 1 : season_num == 2
    season_num == 1 ? season = "Mesotrophic" : season = "Oligotrophic"

    rstar_b, _, _ = rstar_analysis(fsaven, season_num)

    if ds["pulse"][:] == 1
        N, P, Z, B, D, O = get_endpoints(["n", "p", "z", "b", "d", "o"], ds)
    else
        N, P, Z, B, D, O = mean_over_time(["n", "p", "z", "b", "d", "o"], ds, season_num)
    end
    
    f1, dir1 = plot_stacked_dar(N, P, Z, B, D, O, zc, filename)
    println("saving to: $(dir1)/$(filename)_dar.png")
    savefig(f1,"$(dir1)/$(filename)_dar.png")

    f2, dir2 = plot_individual_dar(P, Z, B, D, N, zc, filename)
    println("saving to: $(dir2)/$(filename)_dar.png")
    savefig(f2,"$(dir2)/$(filename)_dar.png")

    f3, dir3 = plot_combined(P, Z, B, D, N, zc, filename, rstar_b)
    println("saving to: $(dir3)/$(filename).png")
    savefig(f3,"$(dir3)/$(filename).png")

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



function plot_individual_dar(P, Z, B, D, N, zc, filename)

    parent_folder = "results/plots/bmass_individual/"
    dir = check_subfolder_exists(filename, parent_folder)

    P, Z, B = set_extinct_to_zero(P), set_extinct_to_zero(Z), set_extinct_to_zero(B)

    yl=(-890.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls=7
    ab=0.7
    lfs = 8
    ls2=5

    dcols = ["teal", "lightblue1", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "red2"]
    bcols = ["teal", "lightblue1", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "coral", "grey", "lime", "orchid", "pink2"]
    tls = ["POM", "POM Consumers", "DOM Lab", "DOM S.Lab", "DOM Mid", "DOM S.Rec", "DOM Rec", "Phyto"]

    p1 = plot(D[:,1], -zc, lc=dcols[1], lw=ls2, linestyle=:dot, grid=false,  label=" Most Labile", xrotation=45, ylimits=yl, title=tls[1], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel="")
    plot!(D[:,2], -zc, lc=dcols[2], lw=ls2, linestyle=:dot, label=" ", labelfontsize=lfs, legend=lg)
    plot!(D[:,3], -zc, lc=dcols[3], lw=ls2, linestyle=:dot, label=" Least Labile", labelfontsize=lfs, legend=lg)

    p2 = plot(B[:,1], -zc, lc=bcols[1], lw=ls, grid=false,  label="", xrotation=45, ylimits=yl, title=tls[2], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel="")
    plot!(B[:,2], -zc, lc=bcols[2], lw=ls, label="", labelfontsize=lfs , legend=lg)
    plot!(B[:,3], -zc, lc=bcols[3], lw=ls, label="", labelfontsize=lfs , legend=lg)

    p3 = plot(D[:,4], -zc, lc=dcols[9], lw=ls2, linestyle=:dot, grid=false, label=" D4", xrotation=45, ylimits=yl, title=tls[3], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel="")
    plot!(B[:,4], -zc, lc=bcols[4], lw=ls, label=" Cop", labelfontsize=lfs , legend=lg)
    plot!(B[:,9], -zc, lc=bcols[9], lw=ls, label=" Oli", labelfontsize=lfs , legend=lg)

    p4 = plot(D[:,5], -zc, lc=dcols[9], lw=ls2, linestyle=:dot, grid=false, label=" D5", xrotation=45, ylimits=yl, title=tls[4], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel="")
    plot!(B[:,5], -zc, lc=bcols[5], lw=ls, label=" Cop", labelfontsize=lfs , legend=lg)
    plot!(B[:,10], -zc, lc=bcols[10], lw=ls, label=" Oli.", labelfontsize=lfs , legend=lg)

    p5 = plot(D[:,6], -zc, lc=dcols[9], lw=ls2, linestyle=:dot, grid=false, label=" D6", xrotation=45, ylimits=yl, title=tls[5], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mmol N/m^3")
    plot!(B[:,6], -zc, lc=bcols[6], lw=ls, label=" Cop", labelfontsize=lfs , legend=lg)
    plot!(B[:,11], -zc, lc=bcols[11], lw=ls, label=" Oli", labelfontsize=lfs , legend=lg)

    #NOTE changed D7 and D8 to remove resource line to see biomasss better
    p6 = plot(B[:,7], -zc, lc=bcols[7], lw=ls, grid=false, label=" Cop", xrotation=45, ylimits=yl, title=tls[6], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")
    plot!(B[:,12], -zc, lc=bcols[12], lw=ls, label=" Oli", labelfontsize=lfs , legend=lg)

    p7 = plot(B[:,8], -zc, lc=bcols[8], lw=ls, grid=false, label=" Cop", xrotation=45, ylimits=yl, title=tls[7], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol N/m^3")
    plot!(B[:,13], -zc, lc=bcols[13], lw=ls, label=" Oli", labelfontsize=lfs , legend=lg)

    zcp = zc[1:40]
    p8 = plot(P[1:40,1], -zcp, lc=ncols[1], grid=false, lw=ls, label=" Oli", xrotation=45, title=tls[8], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, xlabel=L" mmol N/m^3", alpha=ab)
    plot!(P[1:40,2], -zcp, lc=pcols[2], lw=ls, label=" ",labelfontsize=lfs , legend=lg, alpha=ab)    
    plot!(P[1:40,3], -zcp, lc=pcols[3], lw=ls, label=" ",labelfontsize=lfs , legend=lg, alpha=ab)
    plot!(P[1:40,4], -zcp, lc=pcols[4], lw=ls, label=" ",labelfontsize=lfs , legend=lg, alpha=ab)
    plot!(P[1:40,5], -zcp, lc=pcols[5], lw=ls, label=" ",labelfontsize=lfs , legend=lg, alpha=ab)
    plot!(P[1:40,6], -zcp, lc=pcols[6], lw=ls, label=" Cop",labelfontsize=lfs , legend=lg, alpha=ab)


    f = plot(p1, p2, p3, p4, p5, p6, p7, p8,
            layout = [1 1 1 1 ; 1 1 1 1],
            fg_legend = :transparent,
            size=(700,600),
            # plot_title = "$season $type",
        )

    
    return f, dir

end

function plot_combined(P, Z, B, D, N, zc, filename, rstar)

    parent_folder = "results/plots/combined/"
    dir = check_subfolder_exists(filename, parent_folder)

    P, Z, B = set_extinct_to_zero(P), set_extinct_to_zero(Z), set_extinct_to_zero(B)

    yl=(-890.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls = 9
    ls2 = 4
    ls3 = 7
    ab = 0.7
    lfs = 7
    ls2 = 4
    xtfs = 8

    dcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "honeydew3"]
    bcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
    rscols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
    tls = ["POM", "POM Consumers", "DOM Lab", "DOM S.Lab", "DOM Mid", "DOM S.Rec", "DOM Rec", "Phyto"]

    # fig1 = Array{Plots.Plot, 1}(undef, 2);
    # fig1[1] =   plot(D[1:89, 1], -zc, lw=ls, lc=dcols[1], label=" Most Labile", legendfontsize=lfs,
    #             ylabel="Depth (m)", xlabel=L"log(mmol N/m^3)", xrotation=45, title =tls[1], titlefontsize=tfs, grid=false, border=:box, legend=lg, 
    #             xtickfontsize=xtfs, xscale=:log10, alpha=ab)
    #             plot!(D[1:89, 2], -zc, lw=ls, lc=dcols[2], label=" ", legendfontsize=lfs, alpha=ab)
    #             plot!(D[1:89, 3], -zc, lw=ls, lc=dcols[3], label=" Least Labile", legendfontsize=lfs, alpha=ab)
    #             plot!(rstar[1][1:89], -zc, lw=ls2, label=" R*B1", linestyle=:dot, lc=rscols[1])
    #             plot!(rstar[2][1:89], -zc, lw=ls2, label=" R*B2", linestyle=:dot, lc=rscols[2])
    #             plot!(rstar[3][1:89], -zc, lw=ls2, label=" R*B3", linestyle=:dot, lc=rscols[3])

    # fig1[2] =   plot(B[1:89, 1], -zc, lw=ls, lc=bcols[1], label=" B1", legendfontsize=lfs,
    #             ylabel="", xlabel=L" mmol N/m^3", xrotation=45, title=tls[2], titlefontsize=tfs, grid=false, border=:box, legend=lg,
    #             yformatter=Returns(""), xtickfontsize=xtfs)               
    #             plot!(B[1:89, 2], -zc, lw=ls, lc=bcols[2], label=" B2", legendfontsize=lfs)
    #             plot!(B[1:89, 3], -zc, lw=ls, lc=bcols[3], label=" B3", legendfontsize=lfs)

    
    # f1 = plot(fig1..., 
    # fg_legend = :transparent,
    # layout = (1,2),
    # # size=(800,600),
    # )

    fig1 = Array{Plots.Plot, 1}(undef, 12);
    fig1[1] =   plot(D[1:89, 4], -zc, lw=ls3, lc=dcols[9], label=" DOM1", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L"log(mmol N/m^3)", 
                xrotation=45, title=tls[3], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[4][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[4])
                plot!(rstar[9][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[9])

    fig1[2] =   plot(D[1:89, 5], -zc, lw=ls3, lc=dcols[9], label=" DOM2", legendfontsize=lfs, yformatter=Returns(""), xlabel=L"log(mmol N/m^3)", 
                xrotation=45, title=tls[4], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[5][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[5])
                plot!(rstar[10][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[10])

    fig1[3] =   plot(D[1:89, 6], -zc, lw=ls3, lc=dcols[9], label=" DOM3", legendfontsize=lfs, yformatter=Returns(""), xlabel=L"log(mmol N/m^3)", 
                xrotation=45, title=tls[5], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[6][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[6])
                plot!(rstar[11][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[11])

    fig1[4] =   plot(D[1:89, 7], -zc, lw=ls3, lc=dcols[9], label=" DOM4", legendfontsize=lfs,
                yformatter=Returns(""), xlabel=L"log(mmol N/m^3)", xrotation=45, title=tls[6], titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[7][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[7])
                plot!(rstar[12][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[12])

    fig1[5] =   plot(D[1:89, 8], -zc, lw=ls3, lc=dcols[9], label=" DOM5", legendfontsize=lfs,
                yformatter=Returns(""), xlabel=L"log(mmol N/m^3)", xrotation=45, title=tls[7], titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[8][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[8])
                plot!(rstar[13][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[13])

    fig1[6] =   plot(D[1:89, 1], -zc, lw=ls, lc=dcols[1], label=" Lab", legendfontsize=lfs,
                yformatter=Returns(""), xlabel=L"log(mmol N/m^3)", xrotation=45, title =tls[1], titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(D[1:89, 2], -zc, lw=ls, lc=dcols[2], label=" S.Lab", legendfontsize=lfs, alpha=ab)
                plot!(D[1:89, 3], -zc, lw=ls, lc=dcols[3], label=" Rec", legendfontsize=lfs, alpha=ab)
                plot!(rstar[1][1:89], -zc, lw=ls2, label=" R*B1", linestyle=:dot, lc=rscols[1])
                plot!(rstar[2][1:89], -zc, lw=ls2, label=" R*B2", linestyle=:dot, lc=rscols[2])
                plot!(rstar[3][1:89], -zc, lw=ls2, label=" R*B3", linestyle=:dot, lc=rscols[3])


    fig1[7] =   plot(B[1:89, 4], -zc, lw=ls, lc=bcols[4], label=" Cop", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L" mmol N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 9], -zc, lw=ls, lc=bcols[9], label=" Oli", legendfontsize=lfs)

    fig1[8] =   plot(B[1:89, 5], -zc, lw=ls, lc=bcols[5], label=" Cop", legendfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 10], -zc, lw=ls, lc=bcols[10], label=" Oli", legendfontsize=lfs)

    fig1[9] =   plot(B[1:89, 6], -zc, lw=ls, lc=bcols[6], label=" Cop", legendfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 11], -zc, lw=ls, lc=bcols[11], label=" Oli", legendfontsize=lfs)

    fig1[10] =   plot(B[1:89, 7], -zc, lw=ls, lc=bcols[7], label=" Cop", legendfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 12], -zc, lw=ls, lc=bcols[12], label=" Oli", legendfontsize=lfs)

    fig1[11] =   plot(B[1:89, 8], -zc, lw=ls, lc=bcols[8], label=" Cop", legendfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 13], -zc, lw=ls, lc=bcols[13], label=" Oli", legendfontsize=lfs)

    fig1[12] =  plot(B[1:89, 1], -zc, lw=ls, lc=bcols[1], label=" B1", legendfontsize=lfs,
                ylabel="", xlabel=L" mmol N/m^3", xrotation=45, titlefontsize=tfs, grid=false, border=:box, legend=lg,
                yformatter=Returns(""), xtickfontsize=xtfs)               
                plot!(B[1:89, 2], -zc, lw=ls, lc=bcols[2], label=" B2", legendfontsize=lfs)
                plot!(B[1:89, 3], -zc, lw=ls, lc=bcols[3], label=" B3", legendfontsize=lfs)


    f1 = plot(fig1..., 
    fg_legend = :transparent,
    layout = (2,6),
    size=(900,700),
    )

    return f1, dir

end


# fsaven = "results/outfiles/Wi50y_231202_15:33_6P3Z13B8D.nc"
# fsaven = "results/outfiles/Wi100y_231202_16:19_6P3Z13B8D.nc"
# fsaven = "results/outfiles/Wi50y_231202_16:48_6P3Z13B8D.nc"
# fsaven = "results/outfiles/Wi100y_231202_17:10_6P3Z13B8D.nc"

# fsaven = "results/outfiles/Wi50y_231202_23:38_6P3Z13B8D.nc"
# fsaven="results/outfiles/Wi100y_231203_10:39_6P3Z13B8D.nc" #winter steady
# fsaven="results/outfiles/Wi100y_231203_11:03_6P3Z13B8D.nc" #winter pulse
# fsaven="results/outfiles/Su100y_231203_14:48_6P3Z13B8D.nc" #summer steady
# fsaven="results/outfiles/Su100y_231203_19:58_6P3Z13B8D.nc"  #summer pulse

# fsaven = "results/outfiles/Wi50y_231205_21:08_6P3Z13B8D.nc" # winter, steady, 10,10,10 sinking
# fsaven="results/outfiles/Wi50y_231205_22:21_6P3Z13B8D.nc" # winter, 10,10,10, mlz=40
# fsaven="results/outfiles/Wi50y_231205_23:29_6P3Z13B8D.nc" # winter, 10,10,10, mlz=80

# fsaven="results/outfiles/Wi50y_240108_16:48_6P3Z13B8D.nc"
# plot_state_vars_dar(fsaven, 1)
# fsaven="results/outfiles/Su50y_240104_14:04_6P3Z13B8D.nc"
# plot_state_vars_dar(fsaven, 2)