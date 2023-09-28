
using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames
using SparseArrays, LinearAlgebra

#------------------------------------------------------------------------------
#                            PLOT RSTAR B
#------------------------------------------------------------------------------

function plot_rstar(rstar_b, rstar_p, fsaven)

    ds = NCDataset(fsaven)
    N, P, Z, B, D = get_endpoints(["n", "p", "z", "b", "d"], ds)
    BD_competitors = group_BD_competitors(sparse(ds["CM"][:]), get_size([D])[1])

    f = plot_rstar_b(fsaven, rstar_b, B, D, BD_competitors, ds)


end


function group_BD_competitors(Cs, nd)

    competitors = Any[]
    for (i, row) in enumerate(eachrow(Cs))
        for (j, col) in enumerate(eachcol(Cs))
            if Cs[i, j] > 0
                push!(competitors, [i, j])
            end
        end
    end

    return get_competitor_dict(competitors, nd) 

end


function get_competitor_dict(competitors, nd)

    out = Dict(name => Any[] for name in collect(1:1:nd))
    for i in competitors
        for j in keys(out) 
            if i[1] == j 
                push!(out[j], i[2]) 
            end 
        end 
    end

    return out

end

#d = Dict(name => Any[] for name in collect(1:1:nd))

# using IterTools
# for i in groupby(x -> x[1], comp)
#     push!(out, i)
# end

function plot_rstar_b(fsaven, rstar, B, D, competitors, ds)

    nb = get_size([B])[1]
    nd = get_size([D])[1]

    H = 500
    dz = 10
    zc = [dz/2:dz:(H-dz/2)]
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tl=:right
    tfs=22

    p1 = plot(rstar[1][1:50], -zc, lw=ls, lc=bcols[1], label=" B1", ylabel="Depth (m)", xrotation=45, xguidefontsize=12, 
    xlabel="", border=:box, legend=lg, xscale=:log10, title="D1", title_loc=tl, titlefontsize=tfs, legendfontsize=lfs)
    plot!(Dw[1:50, 1], -zc, lw=ls, lc=dcols[1], linestyle=:dot,label=" D1", alpha=ab, legendfontsize=lfs)

    p2 = plot(bw[1:50, 1], -zc, lw=ls, lc=bcols[1], label=" B1", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xscale=:log10, legendfontsize=lfs)

    # p3 = plot(rstar[2][1:50], -zc, lw=ls, lc=bcols[2], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    # border=:box, legend=lg, xscale=:log10, title="D2", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    p3 = plot(rstar[6][1:50], -zc, lw=ls, lc=bcols[6], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D2", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    plot!(rstar[10][1:50], -zc, lw=ls, lc=bcols[10],  label="")
    plot!(Dw[1:50, 2], -zc, lw=ls, lc=dcols[2], linestyle=:dot, alpha=ab, label=" D2", legendfontsize=lfs)

    p4 = plot(bw[1:50, 2], -zc, lw=ls, lc=bcols[2], label=" B2", linestyle=:dash, xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab_ext, legendfontsize=lfs)
    plot!(bw[1:50, 6], -zc, lw=ls, lc=bcols[6],  label=" B6", alpha=ab, legendfontsize=lfs)
    plot!(bw[1:50, 10], -zc, lw=ls, lc=bcols[10],  label=" B10", alpha=ab, legendfontsize=lfs)

    p5 = plot(rstar[3][1:50], -zc, lw=ls, lc=bcols[3], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D3", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    # plot!(rstar[7][1:50], -zc, lw=ls, lc=bcols[7], label="")
    plot!(rstar[11][1:50], -zc, lw=ls, lc=bcols[11], label="")
    plot!(Dw[1:50, 3], -zc, lw=ls, lc=dcols[3], linestyle=:dot, alpha=ab, label=" D3", legendfontsize=lfs)

    p6 = plot(bw[1:50, 3], -zc, lw=ls, lc=bcols[3], label=" B3", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab,legendfontsize=lfs)
    plot!(bw[1:50, 7], -zc, lw=ls, lc=bcols[7], label=" B7", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
    plot!(bw[1:50, 11], -zc, lw=ls, lc=bcols[11], label=" B11", alpha=ab, legendfontsize=lfs)

    p7 = plot(rstar[8][1:50], -zc, lw=ls, lc=bcols[8], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D4", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    # p7 = plot(rstar[4][1:50], -zc, lw=ls, lc="coral4", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    # border=:box, legend=lg, xscale=:log10, title="D4", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.001, 0.02))
    # plot!(rstar[8][1:50], -zc, lw=ls, lc=bcols[8], label="")
    # plot!(rstar[12][1:50], -zc, lw=ls, lc="azure4", label="")
    plot!(Dw[1:50, 4], -zc, lw=ls, lc=dcols[4], linestyle=:dot, alpha=0.6, label=" D4", legendfontsize=lfs)

    p8 = plot(bw[1:50, 4], -zc, lw=ls, lc=bcols[4], linestyle=:dash, label=" B4", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab_ext,  legendfontsize=lfs)
    plot!(bw[1:50, 8], -zc, lw=ls, lc=bcols[8], label=" B8", alpha=ab, legendfontsize=lfs)
    plot!(bw[1:50, 12], -zc, lw=ls, lc=bcols[12], linestyle=:dash, label=" B12", alpha=ab_ext, legendfontsize=lfs)

    p9 = plot(rstar[5][1:50], -zc, lw=ls, lc=bcols[5], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D5", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    plot!(rstar[9][1:50], -zc, lw=ls, lc=bcols[9], label="")
    # plot!(rstar[13][1:50], -zc, lw=ls, lc="turquoise1", label="")
    plot!(Dw[1:50, 5], -zc, lw=ls, lc=dcols[5], linestyle=:dot, alpha=0.6, label=" D5", legendfontsize=lfs)

    p10 = plot(bw[1:50, 5], -zc, lw=ls, lc=bcols[5], label=" B5", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, legendfontsize=lfs)
    plot!(bw[1:50, 9], -zc, lw=ls, lc=bcols[9], label=" B9", alpha=ab, legendfontsize=lfs)
    plot!(bw[1:50, 13], -zc, lw=ls, lc=bcols[13], label=" B13", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
    
    winter = plot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
    fg_legend = :transparent,
    layout = (1,10),
    size=(1600,900),
    plot_xaxis="Depth (m)"
    # xlabel = "R*",
    # plot_title="Winter", 
    # plot_titlefontsize = 20,
    # titlefontsize=tfs, titlelocation=:center, 
    )

    savefig(f,"/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rsB/rsB_$(fsaven).png")
    
    return f

end