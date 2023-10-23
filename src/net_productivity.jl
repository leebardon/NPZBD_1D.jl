using NCDatasets, SparseArrays
# using CairoMakie
using LinearAlgebra, Statistics, CSV, DataFrames

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

    mort = zeros(ngrid, nx, days) 
    for i in range(1, nx)
        mort[:,i,:] += (m_l[i] .+ m_q[i] .* Biomass[:,i,:])
    end

    return dropdims(sum(mort, dims=2), dims=2)

end


function uptake_p(P, N, ds)

    vmax_ij = ds["vmax_ij"][:]
    Kp_ij = ds["Kp_ij"][:]
    # temp_fun = get_temp_mod(ds)
    temp_fun = Matrix(CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/win_temp_mod.csv", DataFrame, header=false))
    light = Matrix(get_light())
    CmP = ds["CMp"][:]
    N = dropdims(N, dims=2)

    # uptake = Array{Float64, 2}(undef, 30, 366) 
    uptake = zeros(length(P[:,1,1]), length(P[1,1,:]))
    for i in 1:length(P[1,:,1])
        uptake += P[:,i,:] .* vmax_ij[i] .* (N ./ (N .+ Kp_ij[i])) 
        # uptake += P[:,i,:] .* temp_fun[1:20] .* vmax_ij[i] .* min.(N ./ (N .+ Kp_ij[i]), light[1:20] ./ (light[1:20] .+ 10.0))
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

    grazing_total = zeros(ngrid, days) 
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

    ngrid = size(mortality, 1)
    days = size(mortality, 2)
    loss = zeros(ngrid, days) 

    loss += mortality .+ grazing

    return loss
end


function get_monthly_means(data)
    # 20 ts per day * 30 days = 600 
    depth = size(data, 1)
    num_months = 12
    monthly_means = zeros(depth, num_months)
    grouped = mean.(Iterators.partition(eachcol(data), 600))

    for i in range(1, num_months)
        monthly_means[:,i] = grouped[i]
    end

    # monthly_means = ifelse.(monthly_means .< 0.01, 0, monthly_means)
    
    return monthly_means 

end


function get_npp(fsaven)

    ds = NCDataset(fsaven)
    P, Z, N = get_final_year(ds, ["p", "z", "n"])

    p, z, n = P[1:20,:,:], Z[1:20,:,:], N[1:20,:,:]
    p_survivors = set_extinct_to_zero(p)
    z_survivors = set_extinct_to_zero(z)
    n_survivors = set_extinct_to_zero(n)

    mort = get_mortality(p_survivors, "P", ds)
    grazing = p_grazing(p_survivors, z_survivors, ds)
    loss = get_loss(mort, grazing, ds)
    uptake = uptake_p(p_survivors, n_survivors, ds)
    
    npp = uptake .- loss
    monthly_mean_npp = get_monthly_means(npp)
    display(monthly_mean_npp)

    npp_heatmaps(fsaven, monthly_mean_npp)

end



fsaven = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_231019_19:45_10P3Z18B8D.nc"


get_npp(fsaven)
# npp = get_npp(fsaven)
# npp_heatmaps(fsaven, npp)

# ds = NCDataset(fsaven)
# P, Z, N = get_final_year(ds, ["p", "z", "n"])

# p, z, n = P[1:20,:,:], Z[1:20,:,:], N[1:20,:,:]
# p_survivors = set_extinct_to_zero(p)
# z_survivors = set_extinct_to_zero(z)
# n_survivors = set_extinct_to_zero(n)

# mort = get_mortality(p_survivors, "P", ds)
# grazing = p_grazing(p_survivors, z_survivors, ds)
# loss = get_loss(mort, grazing, ds)
# uptake = uptake_p(p_survivors, n_survivors, ds)
# npp = uptake .- loss
# monthly_mean_npp = get_monthly_means(npp)