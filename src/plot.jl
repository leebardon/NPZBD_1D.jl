using NCDatasets
using Plots, ColorSchemes
using DataFrames

default(show = true)

# outdir = "/home/lee/Dropbox/Development/NPZBD_0D/"
# outfile = "out_0D_20230611.nc"


function get_size(arr)

    out = Vector{Int}()
    
    for a in arr
        append!(out, size(a,1))
    end

    return out[1], out[2], out[3], out[4], out[5]

end


function plot_results(outdir, outfile)

    ds = NCDataset(outdir*outfile)
    file = replace(outfile, "out_0D_" => "", ".nc" => "")
    n, p, z, b, d, timet = get_plot_variables(ds, ["n", "p", "z", "b", "d", "timet"])

    plot_time_total(n, p, z, b, d, timet, file)
    plot_time_individual_pzb(n, p, z, b, d, timet, file)

end


function plot_time_total(n, p, z, b, d, timet, file)

    all_n, all_p, all_z, all_b, all_d = sum_subtypes([n, p, z, b, d])

    p = plot(timet[:], [all_p[:], all_b[:], all_z[:]], lw=2, lc=[:limegreen :skyblue3 :coral3], label=["P" "B" "Z"], grid=false)
    plot!(timet[:], [all_n[:], all_d[:]], lw=2, linecolor=[:darkgrey :olive], ls=:dash, label=["N" "D"])

    title!("Total NPZBD over time")
    xlabel!("Time (days)")
    ylabel!("Concentration")

    savefig(p,"results/plots/time/total/$(file)_t.pdf")

end


function plot_time_individual_pzb(n, p, z, b, d, timet, file)

    nn, np, nz, nb, nd = get_size([n,p,z,b,d])

    p = plot(timet[:], p[1, :], lw=2, linecolor="limegreen", grid=false, label="Phy ($np)")

    if np > 1
        cp1 = Plots.palette(:greens, np)
        for i in 2:np
            plot!(timet[:], p[i, :], lw=2, palette=cp1, grid=false, label=" ")
        end
    end 

    cp2 = Plots.palette(:blues, nb)
    for j in 1:nb
        j == 1 ? lab="Bac ($nb)" : lab=" "
        plot!(timet[:], b[j, :], lw=2, palette=cp2, grid=false, label=lab)
    end

    cp3 = Plots.palette(:reds, nz)
    for k in 1:nz
        k == 1 ? lab="Zoo ($nz)" : lab=" "
        plot!(timet[:], z[k, :], lw=2, palette=cp3, grid=false, label=lab)
    end

    title!("PZB over time")
    xlabel!("Time (days)")
    ylabel!("Concentration (mmol/m^3)")

    savefig(p,"results/plots/time/individual/$(file)_i.pdf")

end



function sum_subtypes(all_states)

    out = Vector{Any}()

    for s in all_states
        append!(out, [sum(s[:,:], dims=1)])
    end

    return out[1], out[2], out[3], out[4], out[5]
end


function get_plot_variables(ds, vars)

    out = Vector{Any}()

    for v in vars
        append!(out, [ds["$v"]])
    end

    return out[1], out[2], out[3], out[4], out[5], out[6]

end



# plot_results(outdir, outfile)