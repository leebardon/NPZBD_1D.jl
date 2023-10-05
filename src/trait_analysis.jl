# using NCDatasets
# using Plots, ColorSchemes, LaTeXStrings
# using DataFrames, NCDatasets
# using SparseArrays, LinearAlgebra

# include("utils/utils.jl")
# include("utils/save_utils.jl")


function remove_extinct(biomass)

    ex = 10^-6
    biomass .= ifelse.(biomass .<= ex, 0.0, biomass) 

    return biomass
    
end


function copiotrophy_indexing(fsaven, prms=nothing)

    ds = NCDataset(fsaven)

    Fg_p =  [0.08, 0.15, 0.25, 0.45, 0.68, 0.79, 0.84, 0.91]
    Fg_b =  [0.5, 0.15, 0.22, 0.31, 0.45, 0.53, 0.60, 0.68, 0.73, 0.79, 0.84, 0.88, 0.91]

    # prms.pulse == 1 ? get_endpoint_copio(Fg_p, Fg_b, ds) : get_copio_over_time(Fg_p, Fg_b, ds)
    Bt, Pt = get_final_year(ds, ["b", "p"])
    Bt_copio = get_copio_over_time(Bt, Fg_b, ds)

end 


function get_total_biomass(biomass)

    total = sum(biomass, dims=2)
    # B_adj = get_adjusted_biomass(winners[2], ds, get_size([B])[1], Fg_b)
    # copio_b = calc_copiotrophy_index(B_tot, B_adj)

    return total

end 


function get_adjusted_biomass(biomass, n, Fg)

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


function get_copio_over_time(biomass, Fg, ds)

    n = get_size([biomass])[1]
    ngrid = length(biomass[:,1,1])
    ts_tot = 7320

    winners = remove_extinct(biomass)
    tot = total_per_ts(winners, ngrid, ts_tot)
    adj = adj_total_per_ts(winners, n, ngrid, Fg, ts_tot)
    
    copio = Array{Float64, 2}(undef, ngrid, ts_tot)
    for t in range(1, ts_tot)
        copio[:,t] = B_adj[:,t] ./ B_tot[:,t]
    end

    return copio

end


function get_final_year(ds, vars)

    final_yr = Vector{Any}()

    for v in vars
        append!(final_yr, [ds[v][:, :, end-(7320-1):end]])
    end

    return final_yr

end


function total_per_ts(biomass, ngrid, ts_tot)

    total = Array{Float64, 2}(undef, ngrid, ts_tot)
    for t in range(1, ts_tot)
        total[:,t] = get_total_biomass(biomass)
    end

    return total

end


function adj_total_per_ts(biomass, n, ngrid, Fg, ts_tot)

    adj_total = Array{Float64, 2}(undef, ngrid, ts_tot)
    for t in range(1, ts_tot)
        adj_total[:,t] = get_adjusted_biomass(biomass, n, Fg)
    end

    return adj_total

end


# fsaven = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_230928_20:47_8P6Z13B5D.nc"
# indx = copiotrophy_indexing(fsaven)