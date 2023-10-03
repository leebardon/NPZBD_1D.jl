using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames, NCDatasets
using SparseArrays, LinearAlgebra

include("utils/utils.jl")
include("utils/save_utils.jl")


function copiotrophy_indexing_static(fsaven)

    ds = NCDataset(fsaven)
    P, B = get_endpoints(["p", "b"], ds)

    try
        Fg_p = ds["Fg_p"][:]
        Fg_b = ds["Fg_b"][:]
    catch
        Fg_p =  [0.08, 0.15, 0.25, 0.45, 0.68, 0.79, 0.84, 0.91]
        Fg_b =  [0.5, 0.15, 0.22, 0.31, 0.45, 0.53, 0.60, 0.68, 0.73, 0.79, 0.84, 0.88, 0.91]
    end

    winners = remove_extinct([P, B])
    B_tot = sum(winners[2], dims=2)
    B_adj = get_adjusted_biomass(winners[2], ds, get_size([B])[1], Fg_b)
    copio_b = calc_copiotrophy_index(B_tot, B_adj)

    # P_tot = sum(winners[1], dims=2)
    # P_adj = get_adjusted_biomass(P, ds, get_size([P])[1], Fg_p)
    # copio_p = calc_copiotrophy_index(P_tot, P_adj)

    return copio_b
end 



function remove_extinct(biomasses)

    for b in biomasses

        ex = 10^-6
        b .= ifelse.(b .<= ex, 0.0, b)
        
    end

    return biomasses
    
end



function get_adjusted_biomass(biomass, ds, n, Fg)

    ngrid = length(biomass[:,1])
    adj_biomass = zeros(Float64, ngrid, n) 

    for i in 1:n
        adj_biomass[:, i] = biomass[:, i] .* Fg[i]
    end

    return sum(adj_biomass, dims=2)

end


function calc_copiotrophy_index(tot, adj)

    copio = adj ./ tot

    return copio

end


function copio_over_time(fsaven)

    ds = NCDataset(fsaven)
    B, P = get_final_2_years(ds, ["b", "p"])

    try
        Fg_p = ds["Fg_p"][:]
        Fg_b = ds["Fg_b"][:]
    catch
        Fg_p =  [0.08, 0.15, 0.25, 0.45, 0.68, 0.79, 0.84, 0.91]
        Fg_b =  [0.5, 0.15, 0.22, 0.31, 0.45, 0.53, 0.60, 0.68, 0.73, 0.79, 0.84, 0.88, 0.91]
    end

    winners = remove_extinct([P, B])
    B_tot_weekly = get_weekly_total(B)



end


function get_final_2_years(ds, vars)

    final_2yrs = Vector{Any}()

    for v in vars
        append!(final_2yrs, [ds[v][:, :, end-14640:end]])
    end

    return final_2yrs

end


function get_weekly_total(biomass)



end

fsaven = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_230928_20:47_8P6Z13B5D.nc"
indx = copiotrophy_indexing(fsaven)