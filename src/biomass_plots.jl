
using NCDatasets
using Plots, Colors, LaTeXStrings
using DataFrames


function plot_biomasses(fsaven, season_num)

    ds = NCDataset(fsaven)
    # file = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    file = replace(fsaven, ".nc" => "", "results/outfiles/endpoints/" => "")
    season_num == 1 ? season = "Winter" : season = "Summer"
    ep = get_endpoints(ds, ["n", "p", "z", "b", "d"])
    n, p, z, b, d = ep[1], ep[2], ep[3], ep[4], ep[5]

    H = 890
    dz = 10
    zc = [dz/2:dz:(H-dz/2)]
    
    f1 = plot_stacked(p, b, z, d, n, zc)
    f2 = plot_individual(b, d, z, p, zc, season)

    savefig(f1,"results/plots/bmass_stacked/$(file).png")
    savefig(f2,"results/plots/bmass_individual/$(file).png")

end

function plot_stacked(p, b, z, d, n, zc)

    yl=(-500.0, 0)
    bc, dc, pc, nc, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    
    p1 = plot(sum(p, dims = 2), -zc, lw=ls, lc="darkgreen", grid=false, xrotation=45, label="Total P", 
    ylimits=yl, alpha=ab, labelfontsize=lfs, legend=lg)
    plot!(sum(b, dims = 2), -zc, lw=ls, lc="skyblue", label="Total B", alpha=ab, labelfontsize=lfs, legend=lg) 
    plot!(sum(z, dims = 2), -zc, lw=ls, lc="maroon", label="Total Z", alpha=ab, labelfontsize=lfs, legend=lg) 
    p2 = plot(sum(d, dims = 2), -zc, lw=ls, lc="orange", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg)
    p3 = plot(sum(n, dims = 2), -zc, lw=ls, lc="grey", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg)
    # p4 = plot(sum(o, dims = 2), -zc, lw=ls, lc="pink", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, alpha=ab, labelfontsize=lfs)

    f1 = plot(p1, p2, p3,
        linewidth = 2,
        layout = [1 1 1],
        size=(500,400),
        fg_legend = :transparent,
        title = ["Biomass" "D" "N"]
    )

    return f1
end


function plot_individual(b, d, z, p, zc, season)

    yl=(-500.0, 0)
    sizes = get_size([b, d, z, p])
    nb, nd, nz, np = sizes[1], sizes[2], sizes[3], sizes[4]
    tfs = 9

    bc, dc, pc, nc, ab, ab_ext, ls, lfs, lg = get_plot_vars()

    # Plot OM
    tls = ["  Organic Matter", "Z Biomass", "B Biomass", "P Biomass"]
    p1 = plot(d[:,1], -zc, linecolor=dc[1], lw=ls, grid=false,  label=" D1", xrotation=45, ylimits=yl, title=tls[1], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg)
    if nd > 1
        for i in 2:nd
            plot!(d[:,i], -zc, lc=dc[i], lw=ls, label=" D$i", alpha=ab, labelfontsize=lfs, legend=lg)
        end
    end
  
    # Plot zoo
    p2 = plot(z[:,1], -zc, linecolor="black", grid=false, label=" Z1", lw=ls, xrotation=45, ylimits=yl, title=tls[2], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg)
    cmp2 = cgrad(:brg, nz, categorical = true)
    if nz > 1
        for j in 2:nz
            plot!(z[:,j], -zc, palette=cmp2, lw=ls, label=" Z$j", alpha=ab, labelfontsize=lfs, legend=lg)
        end
    end

    # Plot bac
    p3 = plot(b[:,1], -zc, linecolor=bc[1], grid=false, lw=ls, label=" B1", xrotation=45, ylimits=yl, title=tls[3], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg)
    if nb > 1
        for k in 2:nb
            plot!(b[:,k], -zc, lc=bc[k], lw=ls, label=" B$k", alpha=ab, labelfontsize=lfs , legend=lg)
        end
    end

    # Plot phyto
    p4 = plot(p[:,1], -zc, linecolor=pc[1], grid=false, lw=ls, label=" P1", xrotation=45, ylimits=yl, title=tls[4], 
    titlefontsize=tfs, alpha=ab, labelfontsize=lfs, legend=lg)
    if np > 1
        for l in 2:np
            plot!(p[:,l], -zc, lc=pc[l], lw=ls, label=" P$l", alpha=ab, labelfontsize=lfs, legend=lg)
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

    return out

end


function get_plots_vars()

    bc= ["cyan3", "darkorange", "indigo", "coral4", "lightcyan4", "magenta2", "orange4", "seagreen4",
    "darkkhaki", "purple", "crimson",  "azure4", "turquoise1"]
    dc= ["blue3", "black", "maroon", "navy", "brown4"]
    pc = ["olivedrab3", "darkgreen","red4", "cyan4", "purple", "black", "hotpink2", "wheat2" ]
    nc = ["blue2"]
    ab=0.8
    ab_ext=0.8
    ls=5
    lfs=9
    lg=:bottomright
    
    return bc, dc, pc, nc, ab, ab_ext, ls, lfs, lg

end

# outfile = "results/outfiles/endpoints/Su100y_230916_11:25_8P6Z13B5D_ep.nc"
# plot_biomasses(outfile, 2)

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