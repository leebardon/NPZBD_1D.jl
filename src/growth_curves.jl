using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames
using SparseArrays, LinearAlgebra

include("utils/utils.jl")
include("utils/save_utils.jl")


function growth_curves(fsaven)

    ds = NCDataset(fsaven)
    N, P, Z, B, D = get_endpoints(["n", "p", "z", "b", "d"], ds)
    BD_competitors = group_competitors(sparse(ds["CM"][:]), get_size([D])[1])
    # Rx = collect(0.1:0.22:21)[1:89]
    # Rx = collect(0.1:0.165:15)[1:89]
    Rx = collect(0.1:0.11:10)[1:89]



    growth_b = calc_growth_b(B, Rx, ds)
    plot_monod_b(fsaven, growth_b, B, Rx, D, BD_competitors)

end


function calc_growth_b(B, Rx, ds)

    II, JJ = get_nonzero_axes(ds["CM"][:])
    vmax = ds["vmax_ij"][:]
    Km = ds["Km_ij"][:]
    y = ds["y_ij"][:]

    growth = ones(length(B[:,1]), length(B[1,:]))

    for j = axes(II, 1)
        growth[:,j] = y[II[j],JJ[j]] .* (vmax[II[j],JJ[j]] .*  Rx) ./ (Rx .+ Km[II[j],JJ[j]])
    end

    return growth

end

# function plot_monod_b(growth, B, Dx, D, lbl, f_str, loc, lims, cols, type="B")
function plot_monod_b(fsaven, growth, B, Dx, D, competitors)

    parent_folder = "results/plots/growth/"
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    H = 500
    zc = get_zc(H)
    nb = get_size([B])[1]
    nd = get_size([D])[1]

    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tl=:right
    tfs=12
    lfs=6
    xy_tfs=6

    fig1 = Array{Plots.Plot, 1}(undef, nd);
    for d in 1:nd
        for (k, v) in competitors
            if k == d 
                b_di = v
                if d == 1
                    fig1[1] = plot(Dx, growth[:, 1], lw=ls, lc=bcols[1], label=" B1", alpha=ab, legendfontsize=lfs,
                    ylabel="Growth Rate", xlabel="", xrotation=45, title = "D1", titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                    tickfontsize=xy_tfs)
                else
                    b = b_di[1]
                    fig1[d] = plot(Dx, growth[:, b], lw=ls, lc=bcols[b], label=" B$b", alpha=ab, legendfontsize=lfs,
                    ylabel="", xlabel="", xrotation=45, title="D$d", titlefontsize=tfs, grid=false, border=:box, legend=lg,
                    xtickfontsize=xy_tfs)
                end
                if length(b_di) > 1
                    for b in b_di[2:end]
                        plot!(Dx, growth[:, b], -zc, lw=ls, label=" B$b", lc=bcols[b], alpha=ab, tickfontsize=xy_tfs)
                    end
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
                    xtickfontsize=xy_tfs)
                else
                    fig2[d] = plot(B[1:50, b_di[1]], -zc, lw=ls, lc=bcols[b_di[1]], label=" B$(b_di[1])", alpha=ab, legendfontsize=lfs,
                    ylabel="", xlabel="", title="Biomass", titlefontsize=tfs, xrotation=45, grid=false, border=:box, legend=lg,
                    yformatter=Returns(""), xtickfontsize=xy_tfs)
                end
                if length(b_di) > 1
                    for b in b_di[2:end]
                        plot!(B[1:50, b], -zc, lw=ls, label=" B$b", lc=bcols[b], alpha=ab, tickfontsize=xy_tfs)
                    end
                end
            end
        end
    end
    
    f = plot(fig1..., fig2...,
    fg_legend = :transparent,
    layout = (2,nd),
    )

    savefig(f, "$(dir)/gr_$(filename).png")

end

#         p1 = plot(Rx, growth[1], lw=ls, lc=cols[1], xrotation=45, 
#             xlabel="", ylabel="Growth Rate", border=:box, title="Growth Rate on$(lbl[3])", label="")
#             plot!(Rx, growth[2], lw=ls, lc=cols[2], label="")
#             plot!(Rx, [growth[1], growth[2]],
#                 frame=:box,
#                 alpha=a_in,
#                 grid=false,
#                 tickfontsize=5,
#                 lw=ls, ylabel="", xrotation=45, 
#                 lc=[cols[1] cols[2]],
#                 xlabel="", title="", label="", xlim=lims[1], ylim=lims[2],
#                 inset=bbox(loc[1],loc[2],loc[3],loc[4], :bottom, :right),
#                 subplot=2
#             )

#         p2 = plot(biomass[1][1:dz], -zc, lw=ls, lc=cols[1], label=lbl[1], xrotation=45, xlabel="", ylabel="Depth (m)", title="OM Conc.",frame=:box)
#             plot!(biomass[2][1:dz], -zc, lw=ls, lc=cols[2], label=lbl[2])
#             plot!(R[1][1:dz], -zc, lw=ls, lc=cols[3], ls=:dot, label=lbl[3])

#         p3 = plot(Rx, growth[3], lw=ls, lc=cols[1], ylabel="Growth Rate", xrotation=45, 
#             xlabel=L" mmol/m^3", border=:box, title="", label="")
#             plot!(Rx, growth[4], lw=ls, lc=cols[2], label="")
#             plot!(Rx, [growth[3], growth[4]],
#                 frame=:box,
#                 grid=false,
#                 tickfontsize=5,
#                 alpha=a_in,
#                 lw=ls, ylabel="", xrotation=45, 
#                 lc=[cols[1] cols[2]],
#                 xlabel="", title="", label="", xlim=lims[3], ylim=lims[4],
#                 inset=bbox(loc[5],loc[6],loc[7],loc[8], :bottom, :right),
#                 subplot=2
#             )

#         p4 = plot(biomass[3][1:dz], -zc, lw=ls, lc=cols[1], label=lbl[1], xrotation=45, xlabel=L"mmol/m^3", ylabel="Depth (m)", title="",frame=:box)
#             plot!(biomass[4][1:dz], -zc, lw=ls, lc=cols[2], label=lbl[2])
#             plot!(R[2][1:dz], -zc, lw=ls, lc=cols[3], ls=:dot, label=lbl[3])
    

#         f = plot(p1, p2, p3, p4, 
#         fg_legend = :transparent,
#         size=(600,700),
#         layout = l,
#         )

#     savefig(f, "/home/lee/Dropbox/Development/NPZBD_1D/results/plots/growth/$(_run)/$(f_str).png")

#     return f

# end

fsaven = "results/outfiles/Su100y_230928_19:47_8P6Z13B5D.nc"
growth_curves(fsaven)

# winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230827_13:45_4P3Z7B4D_ep.nc")
# summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Su100y_230827_17:10_4P3Z7B4D_ep.nc")

# # Get endpoints 
# Nw, Pw, Zw, Bw, Dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
# Ns, Ps, Zs, Bs, Ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

# growthB, R = calc_growthB(Bw, Dw, winter, "Win")
# plot_growth_over_D([growthB[:, 3], growthB[:, 6]], R, [Bw[:,3], Bw[:,6]], Dw[:,3], [" B3", " B6"], " D3")