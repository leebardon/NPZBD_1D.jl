using NCDatasets, SparseArrays
using CairoMakie
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


function uptake_p(P, N, grid_depth, ds)

    vmax_ij = ds["vmax_ij"][:]
    Kp_ij = ds["Kp_ij"][:]
    temp_fun = Matrix(CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/win_temp_mod.csv", DataFrame, header=false))
    light = Matrix(get_light())
    CmP = ds["CMp"][:]
    N = dropdims(N, dims=2)

    uptake = zeros(length(P[:,1,1]), length(P[1,1,:]))
    for i in 1:length(P[1,:,1])
        # uptake += P[:,i,:] .* vmax_ij[i] .* (N ./ (N .+ Kp_ij[i])) 
        uptake += P[:,i,:] .* temp_fun[1:grid_depth] .* vmax_ij[i] .* min.(N ./ (N .+ Kp_ij[i]), light[1:grid_depth] ./ (light[1:grid_depth] .+ 10.0))
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


function get_moving_averages(depth_arrs, n)

    moving_averages = Any[]

    for d in depth_arrs
        mov_av = [sum(@view d[i:(i+n-1)])/n for i in 1:(length(d)-(n-1))]
        diff = length(d) - length(mov_av)
        top = ones(diff)*NaN
        push!(moving_averages, vcat(top, mov_av))
    end
    
    return moving_averages
    
end


function make_chunks(depth, n)

    c = length(depth) ÷ n
    return [depth[1+c*i:(i == n-1 ? end : c*i+c)] for i = 0:n-1]

end


function get_daily_means(depth_arrs)

    daily_means = Any[]

    for d in depth_arrs
        chunk_means = Any[]
        chunks = make_chunks(d, 360)
        for c in chunks; push!(chunk_means, mean(c)); end
        push!(daily_means, collect(chunk_means))
    end

    return daily_means

end


function get_npp(fsaven, grid_depth)

    ds = NCDataset(fsaven)
    P, Z, N = get_final_year(ds, ["p", "z", "n"])

    p, z, n = P[1:grid_depth,:,:], Z[1:grid_depth,:,:], N[1:grid_depth,:,:]
    p_survivors = set_extinct_to_zero(p)
    z_survivors = set_extinct_to_zero(z)
    n_survivors = set_extinct_to_zero(n)

    # mort = get_mortality(p_survivors, "P", ds)
    # grazing = p_grazing(p_survivors, z_survivors, ds)
    # loss = get_loss(mort, grazing, ds)
    uptake = uptake_p(p_survivors, n_survivors, grid_depth, ds)
    
    # npp = uptake .- loss
    npp = uptake

    monthly_mean_npp = get_monthly_means(npp)
    moving_averages = get_moving_averages([npp[1, :], npp[5, :], npp[10, :], npp[15, :], npp[20, :]], 20)
    daily_means = get_daily_means([npp[1, :], npp[5, :], npp[10, :], npp[15, :], npp[20, :]])

    npp_heatmaps(fsaven, monthly_mean_npp, npp, moving_averages, daily_means)

end



fsaven = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_231019_17:22_10P3Z18B8D.nc"
get_npp(fsaven, 30)
