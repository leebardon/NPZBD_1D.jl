using NCDatasets
using DataFrames, NCDatasets
using SparseArrays, LinearAlgebra

include("utils/utils.jl")
include("utils/save_utils.jl")
include("plotting/heatmaps.jl")


function mort_p(P, ds)

    ngrid = length(P[:,1])
    mort = zeros(Float64, ngrid, n) 

    m_lp = ones(n) * 1e-1  
    m_qp = ones(n) * 0.1 
        for i in range(1, n)
            mort[:, i] += (m_lp[i] .+ m_qp[i] .* P[:,i])
        end

    return mort

end


function uptake_p(P)
    
    uptake = dropdims(sum(biomass, dims=2), dims=2)

    return uptake

end


function loss_p(P)
    
    uptake = dropdims(sum(biomass, dims=2), dims=2)

    return uptake

end


function p_grazing(P, Z, ds, np, nz)

    GrM = ds["GrM"][:]
    grazing = Any[]
    g_max = 1.0
    K_g = 1.0

    grazing_zi = zeros(Float64, length(P[:,1]), np) 
    for z in range(1, nz)
        if sum(GrM[z,1:np]) > 0 
            prey = GrM[z,1:np]' .* P[:,1:end]
            g = g_max .* prey ./ (prey .+ K_g)
            grazing_zi += (g .* Z[:,z] .* GrM[z,1:np]') ./ prey
            grazing_zi = replace!(grazing_zi, NaN => 0.0)
            push!(grazing, grazing_zi)
        end   
    end

        return sum(grazing)

end


function loss(Biomass, mortality, grazing, n)

    ngrid = length(Biomass[:,1])
    loss = zeros(Float64, ngrid, n)

    for i in range(1, n)
        loss[:,i] = mortality[:,i] .+ grazing[:,i]
    end

    return loss
end


function get_npp(fsaven)

    ds = NCDataset(fsaven)
    P, Z = get_final_year(ds, ["p", "z"])

    mort = mort_p()
    grazing = grazing_p()
    loss = loss()
    uptake = uptake_p()

    npp = uptake .- loss

    npp_heatmaps(fsaven, npp)

end