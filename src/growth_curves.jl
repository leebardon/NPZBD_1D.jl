using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames
using SparseArrays, LinearAlgebra

include("utils/utils.jl")
include("utils/save_utils.jl")


function growth_curves(fsaven, season)

    Rx = collect(0.1:0.22:21)[1:89]

    # competitor_pairs = arr

    growthB = calc_growthB(B, Rx, ds)

end


function calc_growthB(B, Rx, ds)

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


function plot_growth_curves(growth, biomass, Rx, R, lbl, _run, f_str, loc, lims, cols, type="B")

    a_in=0.6
    depth = 400
    dz = Int(depth/10)
    zc = get_zc(depth)
    l = @layout [a{0.5h} b{0.3w} ; 
                c{0.5h} d{0.3w}]
    ls=5

    if type == "P"

        
    
    else
        p1 = plot(Rx, growth[1], lw=ls, lc=cols[1], xrotation=45, 
            xlabel="", ylabel="Growth Rate", border=:box, title="Growth Rate on$(lbl[3])", label="")
            plot!(Rx, growth[2], lw=ls, lc=cols[2], label="")
            plot!(Rx, [growth[1], growth[2]],
                frame=:box,
                alpha=a_in,
                grid=false,
                tickfontsize=5,
                lw=ls, ylabel="", xrotation=45, 
                lc=[cols[1] cols[2]],
                xlabel="", title="", label="", xlim=lims[1], ylim=lims[2],
                inset=bbox(loc[1],loc[2],loc[3],loc[4], :bottom, :right),
                subplot=2
            )

        p2 = plot(biomass[1][1:dz], -zc, lw=ls, lc=cols[1], label=lbl[1], xrotation=45, xlabel="", ylabel="Depth (m)", title="OM Conc.",frame=:box)
            plot!(biomass[2][1:dz], -zc, lw=ls, lc=cols[2], label=lbl[2])
            plot!(R[1][1:dz], -zc, lw=ls, lc=cols[3], ls=:dot, label=lbl[3])

        p3 = plot(Rx, growth[3], lw=ls, lc=cols[1], ylabel="Growth Rate", xrotation=45, 
            xlabel=L" mmol/m^3", border=:box, title="", label="")
            plot!(Rx, growth[4], lw=ls, lc=cols[2], label="")
            plot!(Rx, [growth[3], growth[4]],
                frame=:box,
                grid=false,
                tickfontsize=5,
                alpha=a_in,
                lw=ls, ylabel="", xrotation=45, 
                lc=[cols[1] cols[2]],
                xlabel="", title="", label="", xlim=lims[3], ylim=lims[4],
                inset=bbox(loc[5],loc[6],loc[7],loc[8], :bottom, :right),
                subplot=2
            )

        p4 = plot(biomass[3][1:dz], -zc, lw=ls, lc=cols[1], label=lbl[1], xrotation=45, xlabel=L"mmol/m^3", ylabel="Depth (m)", title="",frame=:box)
            plot!(biomass[4][1:dz], -zc, lw=ls, lc=cols[2], label=lbl[2])
            plot!(R[2][1:dz], -zc, lw=ls, lc=cols[3], ls=:dot, label=lbl[3])
    

        f = plot(p1, p2, p3, p4, 
        fg_legend = :transparent,
        size=(600,700),
        layout = l,
        )

    end

    savefig(f, "/home/lee/Dropbox/Development/NPZBD_1D/results/plots/growth/$(_run)/$(f_str).png")

    return f

end

fsaven = "results/outfiles/Wi100y_230923_17:23_8P6Z13B5D.nc"
growth_curves(fsaven, 1)

winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230827_13:45_4P3Z7B4D_ep.nc")
summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Su100y_230827_17:10_4P3Z7B4D_ep.nc")

# Get endpoints 
Nw, Pw, Zw, Bw, Dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
Ns, Ps, Zs, Bs, Ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

growthB, R = calc_growthB(Bw, Dw, winter, "Win")
plot_growth_over_D([growthB[:, 3], growthB[:, 6]], R, [Bw[:,3], Bw[:,6]], Dw[:,3], [" B3", " B6"], " D3")