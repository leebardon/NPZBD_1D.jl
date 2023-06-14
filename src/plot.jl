using NCDatasets
using Plots, ColorSchemes
using DataFrames

default(show = true)

# outdir = "/home/lee/Dropbox/Development/NPZBD_1D/"
# outfile = "out_1D_20230613.nc"

function depth_plots(outdir, outfile)

    ds = NCDataset(outdir*outfile)
    file = replace(outfile, "out_0D_" => "", ".nc" => "")

    n, p, z, b, d, o = get_depth_variables(ds, ["n", "p", "z", "b", "d", "o"])

    H = ds["H"][:]
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]
    
    f1 = plot_depth_profiles(p, b, z, d, n, zc)
    f2 = plot_stacked_bio_nuts(p, b, z, d, n, o, zc)
    f3 = plot_depth_individual(b, d, zc)

    savefig(f1,"results/plots/depth/$(file)_1.pdf")
    savefig(f2,"results/plots/depth/$(file)_2.pdf")
    savefig(f3,"results/plots/depth/$(file)_3.pdf")

end


function plot_depth_profiles(p, b, z, d, n, zc)

    f1 = plot([sum(p, dims = 2), sum(b, dims = 2), sum(z, dims = 2), sum(d, dims = 2), n], -zc; 
    layout = 5, 
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
    plot!(sum(b, dims = 2), -zc, lc="black", label="Total B") 
    plot!(sum(z, dims = 2), -zc, lc="brown", label="Total Z") 
    p2 = plot(sum(d, dims = 2), -zc, grid=false, xrotation=45, label=false)
    p3 = plot(sum(n, dims = 2), -zc, grid=false, xrotation=45, label=false)
    p4 = plot(sum(o, dims = 2), -zc, grid=false, xrotation=45, label=false)
    f2 = plot(p1, p2, p3, p4,
        linewidth = 2,
        layout = [1 1 1 1],
        fg_legend = :transparent,
        title = ["Biomass" "D" "N" "O2"]
    )

    return f2
end


function plot_depth_individual(b, d, zc)
    
    sizes = get_size([b, d])
    nb, nd = sizes[1], sizes[2]

    p1 = plot(d[:,1], -zc, linecolor = "orange", label=" D Pools", ylabel = "Depth (m)", xlabel = "mmol N/m3")
    cp1 = Plots.palette(:reds, nd)
    for i in 2:nd
        plot!(d[:,i], -zc, palette=cp1, grid=false, label=" ")
    end

    p2 = plot(b[:,1], -zc, linecolor = "purple", label = " B Species", xlabel = "mmol N/m3")
    cp2 = Plots.palette(:blues, nb)
    for j in 2:nb
        plot!(b[:,j], -zc, palette=cp2, grid=false, label=" ")
    end
    
    f3 = plot(p1, p2,
            linewidth = 2,
            layout = 2,
            fg_legend = :transparent,
            title = ["Organic matter" "B Biomass"]
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


function get_depth_variables(ds, vars)

    out = Vector{Any}()

    for v in vars
        append!(out, [ds["$v"][:,:,end]])
    end

    return out[1], out[2], out[3], out[4], out[5], out[6]

end

# depth_plots(outdir, outfile)

# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
#NOTE these 'time' plots are ultimately unnecessary for 1D runs
# function get_time_variables(ds, vars)

#     out = Vector{Any}()

#     for v in vars
#         append!(out, [ds["$v"]])
#     end

#     return out[1], out[2], out[3], out[4], out[5], out[6], out[7]

# end

# function time_plots(outdir, outfile)

#     ds = NCDataset(outdir*outfile)
#     file = replace(outfile, "out_0D_" => "", ".nc" => "")

#     H = ds["H"][:]
#     dz = ds["dz"][:]
#     zc = [dz/2:dz:(H-dz/2)]

#     n, p, z, b, d, o, timet = get_time_variables(ds, ["n", "p", "z", "b", "d", "o", "timet"])

#     plot_time_total(n, p, z, b, d, o, zc, timet, file)
#     plot_time_individual_pzb(n, p, z, b, d, o, zc, timet, file)

# end


# function plot_time_total(n, p, z, b, d, o, zc, timet, file)

#     all_n, all_p, all_z, all_b, all_d, all_o = sum_subtypes([n, p, z, b, d, o])

#     p = plot(timet[:], [all_p[:], all_b[:], all_z[:]], lw=2, lc=[:limegreen :skyblue3 :coral3], label=["P" "B" "Z"], grid=false)
#     plot!(timet[:], [all_n[:], all_d[:], all_o[:]], lw=2, linecolor=[:darkgrey :olive], ls=:dash, label=["N" "D" "O"])

#     title!("Total NPZBD over time")
#     xlabel!("Time (days)")
#     ylabel!("Concentration")

#     savefig(p,"results/plots/time/total/$(file)_t.pdf")

# end


# function plot_time_individual_pzb(n, p, z, b, d, o, zc, timet, file)

#     nn, np, nz, nb, nd = get_size([n,p,z,b,d])

#     p = plot(timet[:], p[1, :], lw=2, linecolor="limegreen", grid=false, label="Phy ($np)")

#     if np > 1
#         cp1 = Plots.palette(:greens, np)
#         for i in 2:np
#             plot!(timet[:], p[i, :], lw=2, palette=cp1, grid=false, label=" ")
#         end
#     end 

#     cp2 = Plots.palette(:blues, nb)
#     for j in 1:nb
#         j == 1 ? lab="Bac ($nb)" : lab=" "
#         plot!(timet[:], b[j, :], lw=2, palette=cp2, grid=false, label=lab)
#     end

#     cp3 = Plots.palette(:reds, nz)
#     for k in 1:nz
#         k == 1 ? lab="Zoo ($nz)" : lab=" "
#         plot!(timet[:], z[k, :], lw=2, palette=cp3, grid=false, label=lab)
#     end

#     title!("PZB over time")
#     xlabel!("Time (days)")
#     ylabel!("Concentration (mmol/m^3)")

#     savefig(p,"results/plots/time/individual/$(file)_i.pdf")

# end


# function sum_subtypes(all_states)

#     out = Vector{Any}()

#     for s in all_states
#         append!(out, [sum(s[:,:], dims=2)])
#     end

#     return out[1], out[2], out[3], out[4], out[5], out[6]
# end