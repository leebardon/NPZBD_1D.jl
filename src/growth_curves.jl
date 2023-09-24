using NCDatasets
using Plots, ColorSchemes
using DataFrames, CSV
using SparseArrays, LinearAlgebra
using LaTeXStrings


function calc_growthB(B, D, ds, season)

    II, JJ = get_nonzero_axes(ds["CM"][:])
    vmax = ds["vmax_ij"][:]
    Km = ds["Km_ij"][:]
    y = ds["y_ij"][:]
    temp_fun = get_temp_mod(season)

    R = collect(0.1:0.22:21)[1:89]
    growth = ones(length(B[:,1]), length(B[1,:]))
    # growth = Any[]

    for j = axes(II, 1)
        # uptake = (vmax[II[j],JJ[j]] .*  D[:,II[j]]) ./ (D[:,II[j]] .+ Km[II[j],JJ[j]])
        growth[:,j] = y[II[j],JJ[j]] .* (vmax[II[j],JJ[j]] .*  R) ./ (R .+ Km[II[j],JJ[j]])
        # yield = y[II[j],JJ[j]]
        # push!(growth, uptake .* yield)
    end

    return growth, R

end

function plot_growth_over_D(growth, R, biomass, D, lb, R_str)

    depth = 400
    dz = Int(depth/10)
    zc = get_zc(depth)
    l = @layout [a{0.9h} b{0.3w}]

    n = size(biomass)
    p1 = plot(R, growth[1], lw=4, lc="olivedrab3", label=lb[1], ylabel="Growth Rate", xrotation=45, 
    xlabel=L" D (mmol/m^3)", border=:box, title="Growth vs. $R_str")
    # cmp = cgrad(:coolwarm, categorical = true)
    plot!(R, growth[2], lw=4, lc="red4", label=lb[2])
    # for i = 2:n
    #     plot!(R, growth[i], lw=4, palette=cmp, label=lb[i])
    # end

    p2 = plot(biomass[1][1:dz], -zc, lw=4, lc="olivedrab3", label=lb[1], xrotation=45, xlabel=L"mmol/m^3", ylabel="Depth (m)", title="OM Conc.")
    plot!(biomass[2][1:dz], -zc, lw=4, lc="red4", label=lb[2])
    # for j = 2:n
        # plot!(biomass[j][1:dz], -zc, lw=4, palette=cmp, label=lb[j])
    # end
    plot!(D[1:dz], -zc, lw=4, lc="purple4", ls=:dot, label="$R_str")

    f = plot(p1, p2,
    fg_legend = :transparent,
    size=(500,350),
    layout = l,
    )

    # savefig(f)

    return f

end

function get_zc(depth)

    zc = [10/2:10:(depth-10/2)]

    return zc

end



function get_endpoints(ds, vars)

    out = Vector{Any}()

    for v in vars
        append!(out, [ds["$v"][:,:,end]])
    end

    return out[1], out[2], out[3], out[4], out[5]

end

function get_nonzero_axes(M)

    Cs = sparse(M)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 

function get_temp_mod(season)
    #fit to SPOT data (approx 20 to 4, approx 16 to 4)
    if season == "Win"
        temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/win_temp_mod.csv", DataFrame)
    else
        temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/sum_temp_mod.csv", DataFrame)
    end

    return Matrix(temp_mod)
end


winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230827_13:45_4P3Z7B4D_ep.nc")
summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Su100y_230827_17:10_4P3Z7B4D_ep.nc")

# Get endpoints 
Nw, Pw, Zw, Bw, Dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
Ns, Ps, Zs, Bs, Ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

growthB, R = calc_growthB(Bw, Dw, winter, "Win")
plot_growth_over_D([growthB[:, 3], growthB[:, 6]], R, [Bw[:,3], Bw[:,6]], Dw[:,3], [" B3", " B6"], " D3")