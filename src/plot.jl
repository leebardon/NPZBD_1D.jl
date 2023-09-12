
using NCDatasets
using Plots, Colors, LaTeXStrings
using DataFrames


function depth_plots(fsaven, season_num, years, run_type)

    ds = NCDataset(fsaven)
    file = replace(fsaven, "out_$(years)y_" => "", ".nc" => "", "results/outfiles/" => "")

    season_num == 1 ? season = "Winter" : season = "Summer"

    n, p, z, b, d, o = get_endpoints(ds, ["n", "p", "z", "b", "d", "o"])

    H = ds["H"][:]
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]
    
    f1 = plot_stacked(p, b, z, d, n, o, zc, season)
    f2 = plot_individual(b, d, z, p, zc, run_type, season)

    combined = plot(f1, f2, 
        layout = (2,1),
        size=(700,1000),
        plot_title = season,
    )

    savefig(f1,"results/plots/stacked/$(file)_$(season)_$(years)y.png")
    savefig(f2,"results/plots/individual/$(file)_$(season)_$(years)y.png")
    savefig(combined,"results/plots/combined/$(file)_$(season)_$(years)y.png")

end

function plot_stacked(p, b, z, d, n, o, zc, season)

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


function plot_individual(b, d, z, p, zc, run_type, season)

    yl=(-500.0, 0)
    sizes = get_size([b, d, z, p])
    nb, nd, nz, np = sizes[1], sizes[2], sizes[3], sizes[4]
    tfs = 9

    # Plot OM
    tls = ["  Organic Matter", "Z Biomass", "B Biomass", "P Biomass"]
    if run_type == 3
        p1 = plot(d[:,1], -zc, linecolor="yellow2", lw=3, alpha=0.5, label=L" POM", xrotation=45, ylimits=yl, title=tls[1], titlefontsize=tfs)
        for i in 2:nd
            if i == 2
                plot!(d[:,i], -zc, grid=false, linecolor="goldenrod", lw=3, alpha=0.5, label=L" LL")
            elseif i == 3
                plot!(d[:,i], -zc, grid=false, linecolor="cyan", lw=3, alpha=0.5, label=L" ML")
            else
                plot!(d[:,i], -zc, grid=false, linecolor="teal", lw=3, alpha=0.5, label=L" HL")
            end
        end
    else
        p1 = plot(d[:,1], -zc, linecolor="rosybrown2", lw=3, alpha=0.5, label=L" D1", xrotation=45, ylimits=yl, title=tls[1], titlefontsize=tfs)
        c1, c2 = colorant"tan2", colorant"brown4"
        cmp1 = range(c1, stop=c2, length=nd-1)
        for i in 2:nd
            plot!(d[:,i], -zc, palette=cmp1, grid=false, lw=3, alpha=0.5, label=L" D$i")
        end
    end

    # Plot zoo
    if run_type == 3
        p2 = plot(z[:,1], -zc, linecolor="dodgerblue2", grid=false, label=L" PHY_g", lw=3, alpha=0.5, xrotation=45, ylimits=yl, title=tls[2], titlefontsize=tfs)
        plot!(z[:,2], -zc, linecolor="seagreen", lw=3,alpha=0.5,  label=L" POM_g")
        plot!(z[:,3], -zc, linecolor="red4", lw=3, alpha=0.5, label=L" BAC_g")
    else
        p2 = plot(z[:,1], -zc, linecolor="salmon1", grid=false, label=L" Z1", lw=3, alpha=0.5,  xrotation=45, ylimits=yl, title=tls[2], titlefontsize=tfs)
        c1, c2 = colorant"salmon2", colorant"red4"
        cmp2 = range(c1, stop=c2, length=nz-1)
        for j in 2:nz
            plot!(z[:,j], -zc, palette=cmp2, lw=3, alpha=0.5, label=L" Z$i")
        end
    end

    # Plot bac
    if run_type == 3
        p3 = plot(b[:,1], -zc, linecolor="seagreen", grid=false, lw=3, alpha=0.5, label=L" B1_{POM}", xrotation=45, 
        ylimits=yl, title=tls[3], titlefontsize=tfs)
        for k in 2:nb
            if k == 2 
                plot!(b[:,k], -zc, linecolor="red4", lw=3, alpha=0.5, label=L" B2")
            elseif k == 3
                plot!(b[:,k], -zc, linecolor="darkorange", lw=3, alpha=0.5, label=L" B3")
            elseif k == 4
                plot!(b[:,k], -zc, linecolor="gold4", lw=3, alpha=0.5, label=L" B4")
            elseif k == 5
                plot!(b[:,k], -zc, linecolor="lightcyan4", lw=3, alpha=0.5, label=L" B5")
            elseif k == 6
                plot!(b[:,k], -zc, linecolor="magenta2", lw=3, alpha=0.5, label=L" B6")
            else
                plot!(b[:,k], -zc, linecolor="grey6", lw=3,alpha=0.5,  label=L" B7")
            end
        end
    else
        p3 = plot(b[:,1], -zc, linecolor="plum", grid=false, lw=3,alpha=0.5,  ls=:dot, label="", xrotation=45, ylimits=yl, title=tls[3], titlefontsize=tfs)
        c1, c2 = colorant"plum1", colorant"purple4"
        cmp3 = range(c1, stop=c2, length=nb-1)
        for k in 2:nb
            plot!(b[:,k], -zc, palette=cmp3, lw=3, alpha=0.5, label="")
        end
    end

    # Plot phyto
    p4 = plot(p[:,1], -zc, linecolor="powderblue", grid=false, lw=3, alpha=0.5, label=L" P1", xrotation=45, ylimits=yl, title=tls[4], titlefontsize=tfs)
    # plot!(p[:,2], -zc, linecolor="maroon", lw=3, alpha=0.5, label=L" P2")
    # plot!(p[:,2], -zc, linecolor="hotpink2", lw=3, alpha=0.5, label=L" P3")
    # plot!(p[:,2], -zc, linecolor="rosybrown2", lw=3, alpha=0.5, label=L" P4")
    
    # c1, c2 = colorant"pink", colorant"maroon"
    # cmp4 = range(c1, stop=c2, length=np-2)
    # for l in 2:np
    #     if l == np
    #         plot!(p[:,l], -zc, grid=false, linecolor="black", lw=3, label=" Rate opt")
    #     else
    #         plot!(p[:,l], -zc, palette=cmp4, grid=false, lw=3, label="")
    #     end
    # end

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