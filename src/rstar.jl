using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames, CSV
using SparseArrays, LinearAlgebra, Statistics

include("utils/utils.jl")
include("utils/save_utils.jl")
include("plotting/rstar_plots.jl")


function rstar_analysis(fsaven, season_num)

    global fsaven
    ds = NCDataset(fsaven)
    season_num == 1 ? season = "Winter" : season = "Summer"

    if ds["pulse"][:] == 1
        N, P, Z, B, D = get_endpoints(["n", "p", "z", "b", "d"], ds)
    else
        N, P, Z, B, D = mean_over_time(["n", "p", "z", "b", "d"], ds, season)
    end

    rstar_b, rstar_p, rstar_z = get_rstar(B, P, Z, ds)
    # plot_rstar(rstar_b, rstar_p, rstar_z, fsaven)
    plot_rstar_dar(rstar_b, rstar_p, rstar_z, fsaven)

end 


#-----------------------------------------------------------------------------------
#                                     RSTAR B 
#-----------------------------------------------------------------------------------
function get_rstar(B, P, Z, ds)

    nb, nz, np = get_size([B, Z, P])
    
    mort_b = mortality(B, ds, nb, "B")
    grz_b = b_grazing(B, Z, np, nb, nz, ds)
    loss_b = loss(B, mort_b, grz_b, nb)
    rstar_b = RstarB(loss_b, ds)
    
    mort_p = mortality(P, ds, np, "P")
    grz_p = p_grazing(P, Z, ds, np, nz)
    loss_p = loss(P, mort_p, grz_p, np)
    rstar_p = RstarP(loss_p, ds, np)

    loss_z = mortality(Z, ds, nz, "Z")
    rstar_z = RstarZ(loss_z, ds, nz)

    return rstar_b, rstar_p, rstar_z

end


function mortality(Biomass, ds, n, group)

    ngrid = length(Biomass[:,1])
    mort = zeros(Float64, ngrid, n) 

    if group == "B"
        for i in range(1, n)
            mort[:, i] += (ds["m_lb"][i] .+ ds["m_qb"][i] .* Biomass[:,i])
        end
    elseif group == "P"
        m_lp = ones(n) * 1e-1  
        m_qp = ones(n) * 0.1 
        for i in range(1, n)
            mort[:, i] += (m_lp[i] .+ m_qp[i] .* Biomass[:,i])
        end
    elseif group == "Z"
        m_lz = ones(n) * 1e-2  
        m_qz = ones(n) * 1.0
        for i in range(1, n)
            mort[:, i] += (m_lz[i] .+ m_qz[i] .* Biomass[:,i])
        end
    end

    return mort

end


function b_grazing(B, Z, np, nb, nz, ds)

    GrM = ds["GrM"][:,:]
    g_max = 1.0
    K_g = 1.0
    ngrid = length(B[:,1])

    grazing = Any[]
    grazing_zi = zeros(Float64, ngrid, nb) 
    for z in range(1, nz)
        if sum(GrM[z,np+1:end]) > 0 
            prey = GrM[z,np+1:end]' .*B[:,1:end]
            g = g_max .* prey ./ (prey .+ K_g)
            grazing_zi += (g .* Z[:,z] .* GrM[z,np+1:end]') ./ prey
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


function RstarB(loss, ds)

    umax_ij = ds["umax_ij"][:,:]
    Km_ij = ds["Km_ij"][:,:]
    yield = ds["y_ij"][:,:]
    temp_mod = get_temp_mod(ds)
    II, JJ = get_nonzero_axes(ds["CM"][:,:])

    RS = Any[]
    for j = axes(II, 1)
        push!(RS, Km_ij[II[j],JJ[j]] .* loss[:, j] ./ (yield[II[j],JJ[j]] .* umax_ij[II[j],JJ[j]] .* temp_mod .- loss[:, j]))
    end

    # RS_out = check_for_negatives(RS)

    return RS

end

#-----------------------------------------------------------------------------------
#                                     RSTAR P 
#-----------------------------------------------------------------------------------
function p_grazing(P, Z, ds, np, nz)

    GrM = ds["GrM"][:,:]
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


function RstarP(loss, ds, np)

    vmax_ij = ds["vmax_ij"][:,:]
    Kp_ij = ds["Kp_ij"][:,:]
    temp_mod = get_temp_mod(ds)

    RS = Any[]
    for i in range(1, np)
        push!(RS, Kp_ij[i] .* loss[:, i] ./ (vmax_ij[i] .* temp_mod .- loss[:, i]))
    end

    RS_out = check_for_negatives(RS)

    return RS_out

end



#-----------------------------------------------------------------------------------
#                                     RSTAR Z
#-----------------------------------------------------------------------------------
function RstarZ(loss, ds, nz)

    gmax = ones(nz)*1.0
    K_g = ones(nz)*1.0
    yield = ones(nz)*0.3
    temp_mod = get_temp_mod(ds)

    RS = Any[]
    for i in range(1, nz)
        push!(RS, K_g[i] .* loss[:, i] ./ (yield[i] * gmax[i] .* temp_mod .- loss[:, i]))
    end

    # RS_out = check_for_negatives(RS)

    return RS

end




# fsaven = "results/outfiles/Wi100y_231119_13:54_6P3Z13B8D.nc" # meso pulse 100yrs
# fsaven = "results/outfiles/Wi100y_231119_15:34_6P3Z13B8D.nc" # meso steady 100yrs
# fsaven = "results/outfiles/Su100y_231119_14:25_6P3Z13B8D.nc" #oligo pulse 100yrs
fsaven = "results/outfiles/Su100y_231119_16:05_6P3Z13B8D.nc" #oligo steady 100yrs
# rstar_analysis(fsaven, 1)
rstar_analysis(fsaven, 2)
