using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames, CSV
using SparseArrays, LinearAlgebra

include("utils/utils.jl")

# winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230827_1345.nc")
# summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230827_1710.nc")

function rstar_analysis(ds, season=nothing)

    N, P, Z, B, D = get_endpoints(["n", "p", "z", "b", "d"], ds1)
    np = get_size([P])[1]

    rstar_b = get_rstar(B, P, Z, ds, season)
    rstar_P = get_rstar(B, P, Z, ds, season)

    return rstar_b, rstar_p

end 



#-----------------------------------------------------------------------------------
#                                     RSTAR B 
#-----------------------------------------------------------------------------------
function get_rstar(B, P, Z, ds, season=nothing)

    nb = get_size([B])[1]
    nz = get_size([Z])[1]
    np = get_size([P])[1]
    
    mort_b, mort_p = mortality(B, ds, nb, "B"), mortality(P, ds, np, "P")
    grz_b, grz_p = b_grazing(B, Z, np, nb, nz, ds), p_grazing(P, Z, ds, np, nz)
    loss_b, loss_p = loss(B, mort_b, grz_b, nb), loss(P, mort_p, grz_p, np)
    rstar_b, rstar_p = RstarB(loss_b, ds, season), RstarP(loss_p, ds, np)

    return rstar_b, rstar_p

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
    end

    return mort

end


function b_grazing(B, Z, np, nb, nz, ds)

    GrM = ds["GrM"][:]
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


function RstarB(loss, ds, season)

    season !== nothing ? temp_mod = get_temp_mod(season) : temp_mod = ds["temp_fun"][:]

    vmax_ij = ds["vmax_ij"][:]
    Km_ij = ds["Km_ij"][:]
    yield = ds["y_ij"][:]
    II, JJ = get_nonzero_axes(ds["CM"][:])
    RS = Any[]

    for j = axes(II, 1)
        push!(RS, Km_ij[II[j],JJ[j]] .* loss[:, j] ./ (yield[II[j],JJ[j]] .* vmax_ij[II[j],JJ[j]] .* temp_mod .- loss[:, j]))
    end

    for i in range(1, length(RS))
        RS[i] = check_for_negatives(RS[i])
    end

    return RS

end


function get_temp_mod(season)
    #fit to SPOT data (approx 20 to 4, approx 16 to 4)
    if season == "Win"
        temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/win_temp_mod.csv", DataFrame)
    else
        temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/sum_temp_mod.csv", DataFrame)
    end

    return Matrix(temp_mod)
end

#-----------------------------------------------------------------------------------
#                                     RSTAR P 
#-----------------------------------------------------------------------------------
# function get_rstar_P(P, Z, ds, np, nz, season)
    
#     mort_p = p_mortality(P, ds, np)
#     grz_p = p_grazing(P, Z, ds, np, nz)
#     loss_p = p_loss(mort_p, grz_p, np)
#     RstarP = RstarP(loss_p, ds, np, season)

#     return RstarP

# end


# function p_mortality(P, ds, np)

#     m_lp = 1e-1
#     m_qp = 1e-1
#     mort_p = Any[]
#     for i in range(1, np)
#         push!(mort_p, (m_lp .+ m_qp .* P[:,i]))
#     end

#     return mort_p

# end


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


# function p_loss(mortality, grazing, np)

#     if np == 8
#         loss = zeros(Float64, 89, np)
#         for i in range(1, np)
#             loss[:, i] = mortality[i] .+ grazing[:, i]
#         end
#     else
#         loss = Any[]
#         for i in range(1, np)
#             push!(loss, mortality[i] .+ grazing[i])
#         end
#     end

#     return loss

# end


function RstarP(loss, ds, np, season=nothing)

    season !== nothing ? temp_mod = get_temp_mod(season) : temp_mod = ds["temp_fun"][:]
    umax_ij = ds["umax_ij"][:]
    Kp_ij = ds["Kp_ij"][:]
    
    RS = Any[]
    for i in range(1, np)
        push!(RS, Kp_ij[i] .* loss[:, i] ./ (umax_ij[i] .* temp_mod .- loss[:, i]))
    end

    for i in range(1, length(RS))
        RS[i] = check_for_negatives(RS[i])
    end

    return RS

end

# function get_endpoints(ds, vars)

#     out = Vector{Any}()

#     for v in vars
#         append!(out, [ds["$v"][:,:,end]])
#     end

#     return out[1], out[2], out[3], out[4], out[5]

# end

# function get_nonzero_axes(M)

#     Cs = sparse(M)
#     (II, JJ, _) = findnz(Cs) 
    
#     return II, JJ

# end 

# nw, pw, zw, bw, dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
# ns, ps, zs, bs, ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])


ds1 = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230923_17:23_8P6Z13B5D_ep.nc")
ds2 = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230905_20:05_2P2Z2B2D_ep.nc")

# rs1_new, rs1_old = rstar_analysis(ds1, ds1)
rsp, rsb = rstar_analysis(ds1)
