
using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames
using SparseArrays, LinearAlgebra

#------------------------------------------------------------------------------
#                            PLOT RSTAR B
#------------------------------------------------------------------------------

function plot_rstar(rstar_b, rstar_p, rstar_z, fsaven)

    ds = NCDataset(fsaven)
    N, P, Z, B, D = get_endpoints(["n", "p", "z", "b", "d"], ds)
    BD_competitors = group_interactions( sparse(ds["CM"][:,:]), get_size([D])[1] )
    Z_prey = group_interactions( sparse(ds["GrM"][:,:]), (get_size([B])[1] + get_size([P])[1]) ) 

    plot_rstar_b(fsaven, rstar_b, B, D, BD_competitors, ds)
    plot_rstar_p(fsaven, rstar_p, P, N, ds)
    plot_rstar_z(fsaven, rstar_z, Z, B, P, Z_prey, ds)

end


function plot_rstar_b(fsaven, rstar, B, D, competitors, ds)

    H = 500
    zc = get_zc(H)
    nb, nd = get_size([B, D])

    parent_folder = "results/plots/rstar/"
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")
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
                    xtickfontsize=xtfs)
                else
                    fig1[d] = plot(D[1:50, d], -zc, lw=ls, lc=dcols[d], linestyle=:dot, label=" D$d", alpha=ab, legendfontsize=lfs,
                    ylabel="", xlabel="", xrotation=45, title="R*", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                    yformatter=Returns(""), xtickfontsize=xtfs)
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
    nn, np = get_size([N, P])

    parent_folder = "results/plots/rstar/"
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")
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



function plot_rstar_z(fsaven, rstar, Z, B, P, Z_prey, ds)

    H = 500
    zc = get_zc(H)
    nz, np, nb = get_size([Z, P, B])
    ngrid = length(P[:,1])

    parent_folder = "results/plots/rstar/"
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tl=:right
    tfs=12
    lfs=6
    xtfs=6
    dom_lc = ["grey", "darkkhaki", "seagreen4"]

    if nz == 6
        fig1 = Array{Plots.Plot, 1}(undef, nz)
        summed_prey = zeros(Float64, ngrid, nz) 
        for z in 1:nz
            for (k, v) in Z_prey
                if k == z 
                    z_pb = v
                    if z == 1
                        summed_prey[:, z] = sum(P[:, z_pb[1]:end], dims=2)
                        fig1[z] = plot(summed_prey[1:50, z], -zc, lw=ls, linestyle=:dot, lc="darkgreen", label=" P$(z_pb[1]):$(z_pb[end])", alpha=ab, legendfontsize=lfs,
                        ylabel="Depth (m)", xlabel="", xrotation=45, title = "R*", titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                        xtickfontsize=xtfs)
                        plot!(rstar[z][1:50], -zc, lw=ls, label=" Z$z", lc=zcols[z], alpha=ab)

                    elseif z == 2
                        summed_prey[:, z] = sum(P[:, z_pb[1]:end], dims=2)
                        fig1[z] = plot(summed_prey[1:50, z], -zc, lw=ls, linestyle=:dot, lc="lime", label=" P$(z_pb[1]):$(z_pb[end])", alpha=ab, legendfontsize=lfs,
                        ylabel="", xlabel="", xrotation=45, title="R*", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                        yformatter=Returns(""), xtickfontsize=xtfs)

                        plot!(rstar[z][1:50], -zc, lw=ls, label=" Z$z", lc=zcols[z], alpha=ab)

                    elseif z == 3
                        summed_prey[:, z] = B[:, 1]
                        fig1[z] = plot(summed_prey[1:50, z], -zc, lw=ls, linestyle=:dot, lc="darkorange", label=" POM", alpha=ab, legendfontsize=lfs,
                        ylabel="", xlabel="", xrotation=45, title="R*", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                        yformatter=Returns(""), xtickfontsize=xtfs)

                        plot!(rstar[z][1:50], -zc, lw=ls, label=" Z$z", lc=zcols[z], alpha=ab)

                    else
                        summed_prey[:, z] = sum(B[:, z_pb[1]:end], dims=2)
                        fig1[z] = plot(summed_prey[1:50, z], -zc, lw=ls, lc=dom_lc[(z-3)], linestyle=:dot, label=" DOM", alpha=ab, legendfontsize=lfs,
                        ylabel="", xlabel="", xrotation=45, title="R*", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                        yformatter=Returns(""), xtickfontsize=xtfs)

                        plot!(rstar[z][1:50], -zc, lw=ls, label=" Z$z", lc=zcols[z], alpha=ab) 
                    end

                end
            end
        end
        
        lay_out = (1, nz)
else
    fig1 = Array{Plots.Plot, 1}(undef, nz)

    fig1[1] = plot(rstar[1][1:50], -zc, lw=ls, lc=zcols[1], label=" Z1", alpha=ab, legendfontsize=lfs,
    ylabel="Depth (m)", xlabel="", xrotation=45, title = "R*", titlefontsize=tfs, grid=false, border=:box, legend=lg, 
    xtickfontsize=xtfs)
    plot!(P[1:50, 1], -zc, lw=ls, linestyle=:dot, label=" P1", lc=pcols[1], alpha=ab)

    b = 2
    for z in 2:nz
        if z <= np
            fig1[z] = plot(rstar[z][1:50], -zc, lw=ls, label=" Z$(z)", lc=zcols[z], alpha=ab, legendfontsize=lfs,
                    ylabel="", xlabel="", xrotation=45, title="R*", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                    yformatter=Returns(""), xtickfontsize=xtfs)
                    plot!(P[1:50, z], -zc, lw=ls, linestyle=:dot, label=" P$(z)", lc=pcols[z], alpha=ab)
        elseif z > np && z <= 10
            fig1[z] = plot(rstar[z][1:50], -zc, lw=ls, label=" Z$(z)", lc=zcols[z], alpha=ab, legendfontsize=lfs,
                    ylabel="", xlabel="", xrotation=45, title="R*", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                    yformatter=Returns(""), xtickfontsize=xtfs)
                    plot!(B[1:50, b], -zc, lw=ls, linestyle=:dot, label=" B$(b)", lc=bcols[b], alpha=ab)
                    b += 1
        elseif z == 11
            fig1[z] = plot(rstar[z][1:50], -zc, lw=ls, label=" Z$(z)", lc=zcols[z], alpha=ab, legendfontsize=lfs,
                    ylabel="Depth (m)", xlabel="", xrotation=45, grid=false, border=:box, legend=lg,
                    xtickfontsize=xtfs)
                    plot!(B[1:50, b], -zc, lw=ls, linestyle=:dot, label=" B$(b)", lc=bcols[b], alpha=ab)
                    b += 1
        else 
            fig1[z] = plot(rstar[z][1:50], -zc, lw=ls, label=" Z$(z)", lc=zcols[z], alpha=ab, legendfontsize=lfs,
                    ylabel="", xlabel="", xrotation=45, grid=false, border=:box, legend=lg,
                    yformatter=Returns(""), xtickfontsize=xtfs)
                    plot!(B[1:50, b], -zc, lw=ls, linestyle=:dot, label=" B$(b)", lc=bcols[b], alpha=ab)
                    b += 1
        end
    end

    lay_out = (2, 10)
    sze = (1200, 400)

end

f = plot(fig1..., 
fg_legend = :transparent,
layout = lay_out,
size=sze,
)

savefig(f, "$(dir)/rsZ_$(filename).png")


end

#---------------------------------------------------------------------------------------------------------------------
function plot_rstar_dar(rstar_b, rstar_p, rstar_z, fsaven)

    ds = NCDataset(fsaven)
    N, P, Z, B, D = get_endpoints(["n", "p", "z", "b", "d"], ds)
    BD_competitors = group_interactions( sparse(ds["CM"][:,:]), get_size([D])[1] )
    Z_prey = group_interactions( sparse(ds["GrM"][:,:]), (get_size([B])[1] + get_size([P])[1]) ) 

    plot_rstar_b_dar(fsaven, rstar_b, B, D, BD_competitors, ds)
    # plot_rstar_p(fsaven, rstar_p, P, N, ds)
    # plot_rstar_z(fsaven, rstar_z, Z, B, P, Z_prey, ds)

end


function plot_rstar_b_dar(fsaven, rstar, B, D, competitors, ds)

    H = 890
    zc = get_zc(H)
    nb, nd = get_size([B, D])

    parent_folder = "results/plots/rstar/"
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls=7
    ab=0.6
    lfs = 8
    ls2=3
    xtfs=8

    dcols = ["teal", "lightblue1", "azure4", "red4", "black"]
    bcols = ["teal", "lightblue1", "azure4", "red4", "pink", "black","grey"]

    tls = ["POM", "POM Consumers", "DOM", "DOM Consumers"]
    lab_om=sum([D[:,4], D[:,5], D[:,6]])
    rs_lab_cop=sum([rstar[4], rstar[5], rstar[6]])
    rs_lab_oli=sum([rstar[9], rstar[10], rstar[11]])

    fig1 = Array{Plots.Plot, 1}(undef, 5);
    fig1[1] =   plot(D[1:89, 1], -zc, lw=ls2, lc=dcols[1], linestyle=:dot, label=" POM1", legendfontsize=lfs,
                    ylabel="Depth (m)", xlabel="", xrotation=45, title = "POM Labile", titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                    xtickfontsize=xtfs)
                plot!(rstar[1][1:89], -zc, lw=ls, label=" R* Copio.", lc=bcols[1])

    fig1[2] =   plot(D[1:89, 2], -zc, lw=ls2, lc=dcols[2], linestyle=:dot, label=" POM2", legendfontsize=lfs,
                    ylabel="", xlabel="", xrotation=45, title="POM Semi", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                    yformatter=Returns(""), xtickfontsize=xtfs)               
                plot!(rstar[2][1:89], -zc, lw=ls, label=" R* Semi", lc=bcols[2])

    fig1[3] =   plot(D[1:89, 3], -zc, lw=ls2, lc=dcols[3], linestyle=:dot, label=" POM3", legendfontsize=lfs,
                ylabel="", xlabel="", xrotation=45, title="POM Recalc.", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                yformatter=Returns(""), xtickfontsize=xtfs)               
                plot!(rstar[3][1:89], -zc, lw=ls, label=" R* Oligo.", lc=bcols[3])

    fig1[4] =   plot(lab_om[1:89, :], -zc, lw=ls2, lc=dcols[4], linestyle=:dot, label=" DOM", legendfontsize=lfs,
                ylabel="", xlabel="", xrotation=45, title="DOM Labile", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                yformatter=Returns(""), xtickfontsize=xtfs)               
                plot!(rs_lab_cop[1:89], -zc, lw=ls, label=" R* Copio.", lc=bcols[4], alpha=ab)
                plot!(rs_lab_oli[1:89], -zc, lw=ls, label=" R* Oligo.", lc=bcols[5], alpha=ab)

    fig1[5] =   plot(D[1:89, 7], -zc, lw=ls2, lc=dcols[5], linestyle=:dot, label=" DOM", legendfontsize=lfs,
                ylabel="", xlabel="", xrotation=45, title="DOM Recalc.", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                yformatter=Returns(""), xtickfontsize=xtfs)               
                plot!(rstar[7][1:89], -zc, lw=ls, label=" R* Copio.", lc=bcols[6])
                plot!(rstar[12][1:89], -zc, lw=ls, label=" R* Oligo.", lc=bcols[7])

    
    f = plot(fig1..., 
    fg_legend = :transparent,
    layout = (1,5),
    size=(700,380),
    )

    savefig(f, "$(dir)/rsBD_$(filename).png")

end

    # fig2 = Array{Plots.Plot, 1}(undef, nd);
    # for d in 1:nd
    #     for (k, v) in competitors
    #         if k == d 
    #             b_di = v
    #             if d == 1
    #                 fig2[d] = plot(B[1:50, b_di[1]], -zc, lw=ls, lc=bcols[b_di[1]], label=" B$(b_di[1])", alpha=ab, legendfontsize=lfs,
    #                 ylabel="Depth (m)", xlabel="", title="Biomass", titlefontsize=tfs, xrotation=45, grid=false, border=:box, legend=lg,
    #                 xtickfontsize=xtfs)
    #             else
    #                 fig2[d] = plot(B[1:50, b_di[1]], -zc, lw=ls, lc=bcols[b_di[1]], label=" B$(b_di[1])", alpha=ab, legendfontsize=lfs,
    #                 ylabel="", xlabel="", title="Biomass", titlefontsize=tfs, xrotation=45, grid=false, border=:box, legend=lg,
    #                 yformatter=Returns(""), xtickfontsize=xtfs)
    #             end
    #             if length(b_di) > 1
    #                 for b in b_di[2:end]
    #                     plot!(B[1:50, b], -zc, lw=ls, label=" B$b", lc=bcols[b], alpha=ab)
    #                 end
    #             end
    #         end
    #     end
    # end