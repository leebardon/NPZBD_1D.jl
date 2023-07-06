
function depth_plots(fsaven, season_num, years)

    ds = NCDataset(fsaven)
    file = replace(fsaven, "out_$(years)y_" => "", ".nc" => "", "results/outfiles/" => "")

    season_num == 1 ? season = "win" : season = "sum"

    n, p, z, b, d, o = get_endpoints(ds, ["n", "p", "z", "b", "d", "o"])

    H = ds["H"][:]
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]
    
    f1 = plot_stacked_bio_nuts(p, b, z, d, n, o, zc)
    f2 = plot_depth_individual(b, d, z, p, zc)

    fig = plot(f1, f2, 
        layout = (2,1),
        size=(700,1000),
        # xlabel = "Conc. (mmol N/m3)",
        # ylabel = "Depth (m)"
    )

    # savefig(f1,"results/plots/individual/$(file)_$(season)_$(years)y.pdf")
    # savefig(f2,"results/plots/total/$(file)_$(season)_$(years)y.pdf")
    # savefig(fig,"results/plots/combined/$(file)_$(season)_$(years)y.pdf")

    savefig(f1,"results/plots/individual/$(file)_$(season)30m_$(years)y2.png")
    savefig(f2,"results/plots/total/$(file)_$(season)30m_$(years)y2.png")
    savefig(fig,"results/plots/combined/$(file)_$(season)30m_$(years)y2.png")

end


function plot_stacked_bio_nuts(p, b, z, d, n, o, zc)

    p1 = plot(sum(p, dims = 2), -zc, lc="green", grid=false, xrotation=45, label="Total P")
    plot!(sum(b, dims = 2), -zc, lc="blue", label="Total B") 
    plot!(sum(z, dims = 2), -zc, lc="black", label="Total Z") 
    p2 = plot(sum(d, dims = 2), -zc, lc="orange", ls=:dot, grid=false, xrotation=45, label="")
    p3 = plot(sum(n, dims = 2), -zc, lc="grey", ls=:dot, grid=false, xrotation=45, label="")
    p4 = plot(sum(o, dims = 2), -zc, lc="pink", ls=:dot, grid=false, xrotation=45, label="")
    f1 = plot(p1, p2, p3, p4,
        linewidth = 2,
        layout = [1 1 1 1],
        size=(700,500),
        fg_legend = :transparent,
        title = ["Biomass" "D" "N" "O2"]
    )

    return f1
end


function plot_depth_individual(b, d, z, p, zc)
    
    sizes = get_size([b, d, z, p])
    nb, nd, nz, np = sizes[1], sizes[2], sizes[3], sizes[4]

    p1 = plot(d[:,1], -zc, linecolor = "grey", lw=1, label="", xrotation=45)
    cp1 = Plots.palette(:lajollaS, nd)
    for i in 2:nd
        if i == 3
            plot!(d[:,i], -zc, palette=cp1, grid=false, linecolor = "cyan", lw=3,ls=:dot, label=" Most labile")
        elseif i == 2
            plot!(d[:,i], -zc, palette=cp1, grid=false, linecolor = "red",lw=3,ls=:dot, label=" Least labile")
        else
            plot!(d[:,i], -zc, linecolor="grey", grid=false, lw=1, label="")
        end
    end

    p2 = plot(z[:,1], -zc, linecolor = "black", label="", lw=1,  xrotation=45)
    cp3 = Plots.palette(:bilbaoS, nb)
    for j in 2:nz
        plot!(z[:,j], -zc, linecolor = "black", grid=false, lw=1, label="")
    end

    p3 = plot(b[:,1], -zc, linecolor = "grey", label="", lw=1, xrotation=45)
    cp2 = Plots.palette(:hawaiiS, nb)
    for k in 2:nb
        if k == 6
            plot!(b[:,k], -zc, grid=false, linecolor = "red", lw=3,ls=:dot, label=" Affinity opt.")
        elseif k == 7
            plot!(b[:,k], -zc, grid=false, linecolor = "cyan", lw=3,ls=:dot, label=" Rate opt.")
        else
            plot!(b[:,k], -zc, linecolor = "grey", grid=false, lw=1,label="")
        end
    end

    p4 = plot(p[:,1], -zc, linecolor = "green", label="", lw=1,xrotation=45)
    cp2 = Plots.palette(:bamakoS, nb)
    for l in 2:np
        if l == 5 || l == 4
            plot!(p[:,l], -zc, grid=false, linecolor = "red", lw=3,ls=:dot, label=" Affinity opt.")
        elseif l == 2 || l == 6
            plot!(p[:,l], -zc, grid=false, linecolor = "cyan", lw=3,ls=:dot, label=" Rate opt.")
        else
            plot!(p[:,l], -zc, linecolor = "grey", grid=false, lw=1,label="")
        end
    end

    f2 = plot(p1, p2, p3, p4,
            layout = [1 1 1 1],
            fg_legend = :transparent,
            size=(700,450),
            title = ["Organic Matter" "Z Biomass" "B Biomass" "P Biomass"]
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


using NCDatasets
using Plots, ColorSchemes
using DataFrames

default(show = true)
outfile = "results/outfiles/out_30y_20230704_0228.nc"
depth_plots(outfile, 1, 30)




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