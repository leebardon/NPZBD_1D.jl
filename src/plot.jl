using NCDatasets
using Plots, ColorSchemes
using DataFrames

default(show = true)

# outfile = "results/outfiles/out_27y_20230628_0956.nc"


function depth_plots(fsaven, pulse, years)

    ds = NCDataset(fsaven)
    file = replace(fsaven, "out_$(years)y_" => "", ".nc" => "", "results/outfiles/" => "")

    if pulse == 0 
        season = "" 
    elseif pulse == 1
        season = "win"
    else 
        season = "sum"
    end 

    n, p, z, b, d, o = get_endpoints(ds, ["n", "p", "z", "b", "d", "o"])

    H = ds["H"][:]
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]
    
    f1 = plot_depth_profiles(p, b, z, d, n, o, zc)
    f2 = plot_stacked_bio_nuts(p, b, z, d, n, o, zc)
    f3 = plot_depth_individual(b, d, z, p, zc)

    savefig(f1,"results/plots/all/$(file)_$(season)_$(years)y.pdf")
    savefig(f2,"results/plots/total/$(file)_$(season)_$(years)y.pdf")
    savefig(f3,"results/plots/biomass/$(file)_$(season)_$(years)y.pdf")

end


function plot_depth_profiles(p, b, z, d, n, o, zc)

    f1 = plot([sum(p, dims = 2), sum(b, dims = 2), sum(z, dims = 2), sum(d, dims = 2), n], -zc; 
    layout = 5, 
    size=(600,800),
    xrotation=45,
    linewidth = 2, 
    legend = false, 
    grid=false,
    title = ["Total P" "Total B" "Total Z" "Total D" "Total N"]
    )

    return f1
end


function plot_stacked_bio_nuts(p, b, z, d, n, o, zc)

    p1 = plot(sum(p, dims = 2), -zc, lc="green", grid=false, xrotation=45, label="Total P")
    plot!(sum(b, dims = 2), -zc, lc="blue", label="Total B") 
    plot!(sum(z, dims = 2), -zc, lc="black", label="Total Z") 
    p2 = plot(sum(d, dims = 2), -zc, lc="orange", ls=:dot, grid=false, xrotation=45, label=false)
    p3 = plot(sum(n, dims = 2), -zc, lc="grey", ls=:dot, grid=false, xrotation=45, label=false)
    p4 = plot(sum(o, dims = 2), -zc, lc="pink", ls=:dot, grid=false, xrotation=45, label=false)
    f2 = plot(p1, p2, p3, p4,
        linewidth = 2,
        layout = [1 1 1 1],
        fg_legend = :transparent,
        title = ["Biomass" "D" "N" "O2"]
    )

    return f2
end


function plot_depth_individual(b, d, z, p, zc)
    
    sizes = get_size([b, d, z, p])
    nb, nd, nz, np = sizes[1], sizes[2], sizes[3], sizes[4]

    p1 = plot(d[:,1], -zc, linecolor = "orange", label="", ylabel = "Depth (m)", xlabel = "mmol N/m3", xrotation=45)
    cp1 = Plots.palette(:turbid, nd)
    for i in 2:nd
        plot!(d[:,i], -zc, palette=cp1, grid=false, label="")
    end

    p2 = plot(z[:,1], -zc, linecolor = "black", label = "", xlabel = "mmol N/m3", xrotation=45)
    cp3 = Plots.palette(:solar, nb)
    for j in 2:nz
        plot!(z[:,j], -zc, palette=cp3, grid=false, label="")
    end

    p3 = plot(b[:,1], -zc, linecolor = "blue", label = "", xlabel = "mmol N/m3", xrotation=45)
    cp2 = Plots.palette(:dense, nb)
    for k in 2:nb
        plot!(b[:,k], -zc, palette=cp2, grid=false, label="")
    end

    p4 = plot(p[:,1], -zc, linecolor = "green", label = "", xlabel = "mmol N/m3", xrotation=45)
    cp2 = Plots.palette(:algae, nb)
    for l in 2:np
        plot!(p[:,l], -zc, palette=cp2, grid=false, label="")
    end

    f3 = plot(p1, p2, p3, p4,
            linewidth = 2,
            layout = 4,
            fg_legend = :transparent,
            size=(400,800),
            title = ["Organic Matter" "Z Biomass" "B Biomass" "P Biomass"]
        )
    
    return f3

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

# depth_plots(outfile, 2, 27)
