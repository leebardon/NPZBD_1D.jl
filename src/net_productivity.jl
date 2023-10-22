using NCDatasets, SparseArrays
# using CairoMakie
using LinearAlgebra, Statistics, CSV

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")
include("plotting/heatmaps.jl")


function get_mortality(Biomass, type, ds)

    ngrid = length(Biomass[:,1,1])
    nx = length(Biomass[1,:,1])
    days = length(Biomass[1,1,:])

    if type == "P"
        m_l, m_q = ones(nx) * 1e-1, ones(nx) * 0.1 
    else
        m_l, m_q = ds["m_lb"][:], ds["m_qb"][:]
    end

    mort = Array{Float64, 3}(undef, ngrid, nx, days) 
    for i in range(1, nx)
        mort[:,i,:] += (m_l[i] .+ m_q[i] .* Biomass[:,i,:])
    end

    return dropdims(sum(mort, dims=2), dims=2)

end


function uptake_p(P, N, ds)

    vmax_ij = ds["vmax_ij"][:]
    Kp_ij = ds["Kp_ij"][:]
    temp_fun = get_temp_mod(ds)
    light = Matrix(get_light())
    CmP = ds["CMp"][:]
    N = dropdims(N, dims=2)

    uptake = Array{Float64, 2}(undef, 30, 366) 
    for i in 1:length(P[1,:,1])
        uptake += P[:,i,:] .* vmax_ij[i] .* min.(N ./ (N .+ Kp_ij[i]), light[1:30] ./ (light[1:30] .+ 10.0))
    end

    return uptake

end


function p_grazing(P, Z, ds)

    ngrid = length(P[:,1,1])
    np = length(P[1,:,1])
    nz = length(Z[1,:,1])
    days = length(P[1,1,:])

    GrM = ds["GrM"][:,:]
    g_max = 1.0
    K_g = 1.0

    grazing_total = Array{Float64, 2}(undef, ngrid, days) 
    for z in range(1, nz)
        if sum(GrM[z,1:np]) > 0 
            prey = GrM[z,1:np]' .* P
            g = g_max .* prey ./ (prey .+ K_g)
            grazing_total += (dropdims(sum(g, dims=2), dims=2) .* Z[:,z,:]) ./ dropdims(sum(prey, dims=2), dims=2)
        end   
    end

        return grazing_total

end


function get_loss(mortality, grazing, ds)

    ngrid = length(mortality[:,1])
    days = size(mortality[1,:])[1]
    loss = Array{Float64, 2}(undef, ngrid, days) 

    loss += mortality .+ grazing

    return loss
end


function get_npp(fsaven)

    ds = NCDataset(fsaven)
    P, Z, N = get_final_year(ds, ["p", "z", "n"])

    p, z, n = P[1:30,:,:], Z[1:30,:,:], N[1:30,:,:]
    data_p, data_z, data_n = p[:,:,1:20:end], z[:,:,1:20:end], n[:,:,1:20:end]
    mort = get_mortality(data_p, "P", ds)
    grazing = p_grazing(data_p, data_z, ds)
    loss = get_loss(mort, grazing, ds)
    uptake = uptake_p(data_p, data_n, ds)
    npp = uptake .- loss

    npp_heatmaps(fsaven, npp)

end



fsaven = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_231019_17:22_10P3Z18B8D.nc"

get_npp(fsaven)
