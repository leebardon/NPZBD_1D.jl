
using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames
using SparseArrays, LinearAlgebra

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")

#------------------------------------------------------------------------------
#                            TEST FOR EQUILIBRIUM
#------------------------------------------------------------------------------

function final2(ds, vars)

    final_2yrs = Vector{Any}()

    for v in vars
        append!(final_2yrs, [ds[v][:, :, end-14640:end]])
    end

    return final_2yrs

end


function final2_yr1(ds)

    final2_yr1 = Vector{Any}()

    for i in 1:length(["p", "z", "b", "d"])
        append!(final2_yr1, [ds[i][:, :, 1:7320]])
    end

    return final2_yr1[1], final2_yr1[2], final2_yr1[3], final2_yr1[4]

end


function final2_yr2(ds)

    final2_yr2 = Vector{Any}()

    for j in 1:length(["p", "z", "b", "d"])
        append!(final2_yr2, [ds[j][:, :, end-7320:end]])
    end

    return final2_yr2[1], final2_yr2[2], final2_yr2[3], final2_yr2[4]

end


function get_pct_change(yr1, yr2)

    year1 = set_extinct_to_zero(yr1[:,:,end])
    year2 = set_extinct_to_zero(yr2[:,:,end])

    pct_chg = ((year2 .- year1)./ year1) * 100
    pct_chg .= ifelse.(isnan.(pct_chg), 0.0, pct_chg)

    return year1, year2, pct_chg

end


function plot_equilib(yr1, yr2, group)

    H = 890
    dz = 10
    zc = [dz/2:dz:(H-dz/2)]
    n = get_size([yr1])
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs=9
    al=0.5

    if group == "P"
        cols = pcols
    elseif group == "B"
        cols = bcols
    elseif group == "Z"
        cols = zcols
    elseif group == "D"
        cols = dcols
    else
        cols = ncols
    end

    year1, year2, pct_chg = get_pct_change(yr1, yr2)

    p1 = plot(year1[:,1,end], -zc, grid=false, lw=ls, lc=cols[1], title="Penultimate Year", label="", xrotation=45,  
    titlefontsize=tfs, xlabel=L" mm/m^3", alpha=al, ylabel="Depth (m)")
    for i in 2:n[1]
        plot!(year1[:,i,end], -zc, lw=ls, lc=cols[i], label="", labelfontsize=lfs, alpha=al)
    end

    p2 = plot(year2[:,1,end], -zc, grid=false, lw=ls, lc=cols[1], title="Final Year", label="", xrotation=45, 
    titlefontsize=tfs, xlabel=L" mm/m^3", alpha=al, yformatter=Returns(""),)
    for i in 2:n[1]
        plot!(year2[:,i,end], -zc, lw=ls, lc=cols[i], alpha=al, label="")
    end

    p3 = plot(pct_chg[:, 1], -zc, grid=false, lw=ls, lc=cols[1], labelfontsize=lfs, titlefontsize=tfs, xlabel=L" \% ", 
    title="% Change", alpha=al, yformatter=Returns(""), label=" $(group)1", legend=lg, xrotation=45)
    for i in 2:n[1]
        plot!(pct_chg[:, i], -zc, lw=ls, lc=cols[i], label=" $(group)$(i)", labelfontsize=lfs, alpha=al)
    end

    f = plot(p1, p2, p3,
    layout = [1 1 1],
    fg_legend = :transparent,
    size=(600,450),
    plot_title = "Equilibrium State for $(group)",
    )

    savefig(f, "test.png")

end


# function rstar_plots(fsaven, season)

#     ds = NCDataset(fsaven)
#     file = replace(fsaven, ".nc" => "", "results/outfiles/endpoints/" => "")
#     season_num == 1 ? season = "Winter" : season = "Summer"
#     ep = get_endpoints(["n", "p", "z", "b", "d"], ds)
#     n, p, z, b, d = ep[1], ep[2], ep[3], ep[4], ep[5]

# end


winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_230923_17:23_8P6Z13B5D.nc")
summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Su100y_230923_19:36_8P6Z13B5D.nc")

final2_w = final2(winter, ["p", "z", "b", "d" ])
final2_s = final2(summer, ["p", "z", "b", "d" ])

pw_yr1, zw_yr1, bw_yr1, dw_yr1 = final2_yr1(final2_w)
pw_yr2, zw_yr2, bw_yr2, dw_yr2 = final2_yr2(final2_w)

plot_equilib(pw_yr1, pw_yr2, "P")