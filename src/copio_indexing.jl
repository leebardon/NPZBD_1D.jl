
using DataFrames, NCDatasets
using SparseArrays, LinearAlgebra

include("utils/utils.jl")
include("utils/save_utils.jl")
include("plotting/heatmaps.jl")



function copio_index_analysis(fsaven, prms=nothing)

    ds = NCDataset(fsaven)
    Fg_p =  ds["Fg_p"][:]
    Fg_b =  ds["Fg_b"][:]

    B, P = get_final_year(ds, ["b", "p"])
    copio_b = get_copio_over_time(B, Fg_b)
    copio_p = get_copio_over_time(P, Fg_p)

    copio_heatmaps(fsaven, copio_b, "B")
    copio_heatmaps(fsaven, copio_p, "P")

end 


function get_total_biomass(biomass)

    total = dropdims(sum(biomass, dims=2), dims=2)

    return total

end 


function calc_copiotrophy_index(tot, adj)

    copio = adj ./ tot

    return copio

end


function get_copio_over_time(biomass, Fg)

    winners = set_extinct_to_zero(biomass)
    tot = total_biomass(winners)
    adj = adj_total_biomass(winners, Fg)
    
    copio = Array{Float64, 2}(undef, size(biomass, 1), size(biomass, 3))
    for t in 1:size(biomass, 3)
        copio[:,t] = adj[:,t] ./ tot[:,t]
    end

    return copio

end


# function get_final_year(ds, vars)

#     final_yr = Vector{Any}()

#     for v in vars
#         if v != "o"
#             append!(final_yr, [ds[v][:, :, end-7319:end]])
#         else
#             append!(final_yr, [ds[v][:, end-7319:end]])
#         end
#     end

#     return final_yr

# end


function adj_total_biomass(biomass, Fg)

    adj_biomass = zeros(Float64, size(biomass, 1), size(biomass, 2), size(biomass, 3)) 

    for i in 1:size(biomass, 2)
        adj_biomass[:,i,:] = biomass[:,i,:] .* Fg[i]
    end

    # return adj_biomass
    return dropdims(sum(adj_biomass, dims=2), dims=2)

end


function total_biomass(biomass)
    
    total = dropdims(sum(biomass, dims=2), dims=2)

    return total

end

# fsaven = "results/outfiles/Wi100y_231011_20:23_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi100y_231011_23:28_8P20Z13B5D.nc"
# fsaven = "results/outfiles/Wi100y_231017_01:23_10P3Z21B9D.nc"
# fsaven = "results/outfiles/Wi100y_231017_12:01_10P3Z21B9D.nc"
# fsaven = "results/outfiles/Wi50y_231017_23:00_10P3Z18B8D.nc"
# fsaven = "results/outfiles/Wi50y_231017_23:32_10P3Z18B8D.nc"
# fsaven = "results/outfiles/Wi100y_231017_14:44_10P3Z21B9D.nc"
# fsaven = "results/outfiles/Wi100y_231019_17:22_10P3Z18B8D.nc"
# fsaven = "results/outfiles/Wi100y_231019_19:45_10P3Z18B8D.nc"
# copio = copio_index_analysis(fsaven)

