
using NCDatasets
using Plots, Colors, LaTeXStrings
using DataFrames


function plot_biomass(fsaven, season_num)

    ds = NCDataset(fsaven)
    file = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    season_num == 1 ? season = "Winter" : season = "Summer"
    ep = get_endpoints(ds, ["n", "p", "z", "b", "d", "o"])
    n, p, z, b, d, o = ep[1], ep[2], ep[3], ep[4], ep[5], ep[6]

    H = ds["H"][:]
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]
    
    f1 = plot_stacked(p, b, z, d, n, o, zc)
    f2 = plot_individual(b, d, z, p, zc, season)

    savefig(f1,"results/plots/bmass_stacked/$(file).png")
    savefig(f2,"results/plots/bmass_individual/$(file).png")

end

function plot_stacked(p, b, z, d, n, o, zc)

    yl=(-400.0, 0)
    p1 = plot(sum(p, dims = 2), -zc, lc="green", grid=false, xrotation=45, label="Total P", ylimits=yl)
    plot!(sum(b, dims = 2), -zc, lc="blue", label="Total B") 
    plot!(sum(z, dims = 2), -zc, lc="black", label="Total Z") 
    p2 = plot(sum(d, dims = 2), -zc, lc="orange", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl)
    p3 = plot(sum(n, dims = 2), -zc, lc="grey", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl)
    p4 = plot(sum(o, dims = 2), -zc, lc="pink", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl)

    f1 = plot(p1, p2, p3, p4,
        linewidth = 2,
        layout = [1 1 1 1],
        size=(700,500),
        fg_legend = :transparent,
        title = ["Biomass" "D" "N" "O2"]
    )

    return f1
end


function plot_individual(b, d, z, p, zc, season)

    yl=(-500.0, 0)
    sizes = get_size([b, d, z, p])
    nb, nd, nz, np = sizes[1], sizes[2], sizes[3], sizes[4]
    tfs = 9

    # Plot OM
    tls = ["  Organic Matter", "Z Biomass", "B Biomass", "P Biomass"]
    p1 = plot(d[:,1], -zc, linecolor="rosybrown2", lw=3, label=" D1", xrotation=45, ylimits=yl, title=tls[1], titlefontsize=tfs)
    cmp1 = cgrad(:coolwarm, nd-1, categorical = true)
    if nd > 1
        for i in 2:nd
            plot!(d[:,i], -zc, palette=cmp1, grid=false, lw=3, label=" D$i")
        end
    end
  
    # Plot zoo
    p2 = plot(z[:,1], -zc, linecolor="salmon1", grid=false, label=" Z1", lw=3, xrotation=45, ylimits=yl, title=tls[2], titlefontsize=tfs)
    cmp2 = cgrad(:brg, nz-1, categorical = true)
    if nz > 1
        for j in 2:nz
            plot!(z[:,j], -zc, palette=cmp2, lw=3, label=" Z$j")
        end
    end

    # Plot bac
    p3 = plot(b[:,1], -zc, linecolor="plum", grid=false, lw=3, label=" B1", xrotation=45, ylimits=yl, title=tls[3], titlefontsize=tfs)
    cmp3 = cgrad(:plasma, nb-1)
    if nb > 1
        for k in 2:nb
            plot!(b[:,k], -zc, palette=cmp3, lw=3, label=" B$k")
        end
    end

    # Plot phyto
    p4 = plot(p[:,1], -zc, linecolor="hotpink2", grid=false, lw=3, label=" P1", xrotation=45, ylimits=yl, title=tls[4], titlefontsize=tfs)
    cmp4 = cgrad(:cool, np-1, categorical = true)
    if np > 1
        for l in 2:np
            plot!(p[:,l], -zc, palette=cmp4, grid=false, lw=3, label=" P$l")
        end
    end


    f2 = plot(p1, p2, p3, p4,
            layout = [1 1 1 1],
            fg_legend = :transparent,
            size=(700,450),
            plot_title = season,
        )
    
    return f2

end

# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
function growth_plots()


end










# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------

function get_size(arr)

    out = Vector{Int}()
    
    for a in arr
        append!(out, size(a, 2))
    end

    return out

end


function get_endpoints(ds, vars)

    out = Vector{Any}()

    for v in vars
        append!(out, [ds["$v"][:,:,end]])
    end

    return out[1], out[2], out[3], out[4], out[5], out[6]

end

outfile = "results/outfiles/Su100y_8P6Z13B5D_230915_17:36.nc"
biomass_plots(outfile, 2, 100, 3)

# default(show = true)
# outfile = "results/outfiles/out_100y_20230905_1716.nc"
# depth_plots(outfile, 1, 100, 3)


# function plot_mlz_temp_test(ds1, ds2, ds3, ds4)

#     p1, z1, b1= get_endpoints(ds1, ["p", "z", "b"])
#     p2, z2, b2= get_endpoints(ds2, ["p", "z", "b"])
#     p3, z3, b3= get_endpoints(ds3, ["p", "z", "b"])
#     p4, z4, b4= get_endpoints(ds4, ["p", "z", "b"])

#     H = ds1["H"][:]
#     dz = ds1["dz"][:]
#     zc = [dz/2:dz:(H-dz/2)]

#     pp1 = plot(sum(p1, dims = 2), -zc, lc="green", grid=false, xrotation=45, label="")
#     plot!(sum(b1, dims = 2), -zc, lc="blue", label="") 
#     plot!(sum(z1, dims = 2), -zc, lc="black", label="") 

#     pp2 = plot(sum(p2, dims = 2), -zc, lc="green", grid=false, xrotation=45, label="")
#     plot!(sum(b2, dims = 2), -zc, lc="blue", label="") 
#     plot!(sum(z2, dims = 2), -zc, lc="black", label="") 

#     p3 = plot(sum(p1, dims = 2), -zc, lc="green", grid=false, xrotation=45, label="")
#     plot!(sum(b1, dims = 2), -zc, lc="blue", label="") 
#     plot!(sum(z1, dims = 2), -zc, lc="black", label="") 

#     p4 = plot(sum(p1, dims = 2), -zc, lc="green", grid=false, xrotation=45, label="Total P")
#     plot!(sum(b1, dims = 2), -zc, lc="blue", label="Total B") 
#     plot!(sum(z1, dims = 2), -zc, lc="black", label="Total Z") 

#     f1 = plot(p1, p2, p3, p4,
#         linewidth = 2,
#         layout = [1 1 1 1],
#         size=(700,500),
#         fg_legend = :transparent,
#     )

#     return f1
# end

# ds1 = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_30y_20230703_2230.nc")
# ds2 = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_30y_20230703_2245.nc")
# ds3 = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_30y_20230703_2304.nc")
# ds4 = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_30y_20230703_2321.nc")