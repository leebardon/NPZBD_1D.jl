
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
    BD_competitors = group_competitors(sparse(ds["CM"][:]), get_size([D])[1])

    plot_rstar_b(fsaven, rstar_b, B, D, BD_competitors, ds)
    plot_rstar_p(fsaven, rstar_p, P, N, ds)

end


function plot_rstar_b(fsaven, rstar, B, D, competitors, ds)

    H = 500
    zc = get_zc(H)
    nb = get_size([B])[1]
    nd = get_size([D])[1]

    parent_folder = "results/plots/rstar/"
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tl=:right
    tfs=12
    lfs=6
    xtfs=6

    #TODO add line to remove rstar lines where B is 0/extinct (use set_extinct_to_zero())
    fig1 = Array{Plots.Plot, 1}(undef, nd);
    for d in 1:nd
        for (k, v) in competitors
            if k == d 
                b_di = v
                if d == 1
                    fig1[d] = plot(D[1:50, d], -zc, lw=ls, lc=dcols[d], linestyle=:dot, label=" D$d", alpha=ab, legendfontsize=lfs,
                    ylabel="Depth (m)", xlabel="", xrotation=45, title = "R*", titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                    xtickfontsize=xtfs, xscale=:log10)
                else
                    fig1[d] = plot(D[1:50, d], -zc, lw=ls, lc=dcols[d], linestyle=:dot, label=" D$d", alpha=ab, legendfontsize=lfs,
                    ylabel="", xlabel="", xrotation=45, title="R*", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                    yformatter=Returns(""), xtickfontsize=xtfs, xscale=:log10)
                end
                for b in b_di[1:end]
                    plot!(rstar[b][1:50], -zc, lw=ls, label=" B$b", lc=bcols[b], alpha=ab)
                end
            end
        end
    end

    fig2 = Array{Plots.Plot, 1}(undef, nd);
    for d in 1:nd
        for (k, v) in competitors
            if k == d 
                b_di = v
                if d == 1
                    fig2[d] = plot(B[1:50, b_di[1]], -zc, lw=ls, lc=bcols[b_di[1]], label=" B$(b_di[1])", alpha=ab, legendfontsize=lfs,
                    ylabel="Depth (m)", xlabel="", title="Biomass", titlefontsize=tfs, xrotation=45, grid=false, border=:box, legend=lg,
                    xtickfontsize=xtfs)
                else
                    fig2[d] = plot(B[1:50, b_di[1]], -zc, lw=ls, lc=bcols[b_di[1]], label=" B$(b_di[1])", alpha=ab, legendfontsize=lfs,
                    ylabel="", xlabel="", title="Biomass", titlefontsize=tfs, xrotation=45, grid=false, border=:box, legend=lg,
                    yformatter=Returns(""), xtickfontsize=xtfs)
                end
                if length(b_di) > 1
                    for b in b_di[2:end]
                        plot!(B[1:50, b], -zc, lw=ls, label=" B$b", lc=bcols[b], alpha=ab)
                    end
                end
            end
        end
    end
    
    f = plot(fig1..., fig2...,
    fg_legend = :transparent,
    layout = (2,nd),
    )

    savefig(f, "$(dir)/rsBD_$(filename).png")

end

function plot_rstar_p(fsaven, rstar, P, N, ds)

    #NOTE needs to be merged with plot_rstar_b if we want to use more N's
    H = 200
    zc = get_zc(H)
    nn = get_size([N])[1]
    np = get_size([P])[1]

    parent_folder = "results/plots/rstar/"
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tl=:right
    tfs=12
    lfs=8
    xtfs=8

    fig1 = Array{Plots.Plot, 1}(undef, nn);
    for n in 1:nn
        if n == 1
            fig1[n] = plot(N[1:20, n], -zc, lw=ls, lc=ncols[n], linestyle=:dot, label=" N$n", alpha=ab, legendfontsize=lfs,
            ylabel="Depth (m)", xlabel="", xrotation=45, title = "R*", titlefontsize=tfs, grid=false, border=:box, legend=lg, 
            xtickfontsize=xtfs, xscale=:log10)
        elseif n > 1
            fig1[n] = plot(N[1:20, d], -zc, lw=ls, lc=ncols[n], linestyle=:dot, label=" N$n", alpha=ab, legendfontsize=lfs,
            ylabel="", xlabel="", xrotation=45, title="R*", titlefontsize=tfs, grid=false, border=:box, legend=lg,
            yformatter=Returns(""), xtickfontsize=xtfs, xscale=:log10)
        else
        end

        for p in 1:np
            plot!(rstar[p][1:20], -zc, lw=ls, label=" P$p", lc=pcols[p], alpha=ab)
        end 
    end

    fig2 = Array{Plots.Plot, 1}(undef, nn);
    for p in 1:np
        if p == 1
            fig2[1] = plot(P[1:20, 1], -zc, lw=ls, lc=pcols[1], label=" P1", alpha=ab, legendfontsize=lfs,
            ylabel="", xlabel="", xrotation=45, title = "Biomass", titlefontsize=tfs, grid=false, border=:box, legend=lg, 
            xtickfontsize=xtfs, yformatter=Returns(""))
        else
            plot!(P[1:20, p], -zc, lw=ls, label=" P$p", lc=pcols[p], alpha=ab)
        end
    end

    f = plot(fig1..., fig2...,
    fg_legend = :transparent,
    layout = (nn,2),
    size=(600,400),
    )

    savefig(f, "$(dir)/rsPN_$(filename).png")

end