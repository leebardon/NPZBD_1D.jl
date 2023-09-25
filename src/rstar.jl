# using NCDatasets
# using Plots, ColorSchemes, LaTeXStrings
# using DataFrames
# using SparseArrays, LinearAlgebra

# winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230827_1345.nc")
# summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230827_1710.nc")


function rstar_analysis(P, Z, B, ds)

    RstarB_ij = get_rstar_B(B, Z, ds, nb, nz, season)
    RstarP = get_rstar_P(P, Z, ds, np, nz, season)

    return RstarB_ij, RstarP

end 



#-----------------------------------------------------------------------------------
#                                     RSTAR B 
#-----------------------------------------------------------------------------------
function get_rstar_B(B, Z, ds, nb, nz, season=NaN)
    
    mort_b = b_mortality(B, ds)
    grz_b = b_grazing(B, Z, ds)
    loss_b = b_loss(mort_b, grz_b)
    RstarB_ij = Rstar(loss_b, ds, season)

    return RstarB_ij

end


function b_mortality(B, ds)

    n = get_size([B])

    mort_b = Any[]
    for i in range(1, n)
        push!(mort_b, (ds["m_lb"][i] .+ ds["m_qb"][i] .* B[:,i]))
    end

    return mort_b

end

# function get_prey(GrM, B)
#     B_dom = B[:,2:end]
#     GrM_dom = GrM[4:end,10:end]

#     for (i, row) in enumerate(eachrow(GrM_dom))
#         for (j, col) in enumerate(eachcol(GrM_dom))
#             if GrM_dom[i, j] > 0
#                 prey = GrM[i, j] .* B[:,j]
#             end
#         end
#     end
# end

function b_grazing(B, Z, ds)
    #NOTE this all needs to be generalized / scaled !! 
    GrM = ds["GrM"][:]
    grazing = Any[]
    g_max = 1.0
    K_g = 1.0

    nb = get_size([B])
    nz = get_size([Z])

    #------- for 1N 8P 6Z 13B 5D
    if nb == 13
        # NOTE length(B[:,1]) = prms.ngrid
        grazing_zi = zeros(Float64, length(B[:,1]), nb) 
        np = 8
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

    #------- for 1N 4P 3Z 7B 4D
    if nb==7
        np=4
        prey_POM = GrM[2,5]' .* B[:,1]
        gb_POM = g_max .* prey_POM ./ (prey_POM .+ K_g)
        grz_POM = (gb_POM .* Z[:,2] .* GrM[2,5]' ./ prey_POM) 
        push!(grazing, grz_POM)

        prey = GrM[3,6:end]' .* B[:,2:end]
        for i in range(1, nb-1)
            gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
            grz_i = (gb_i .* Z[:,3] .* GrM[3,i+(np+1)]' ./ prey[:,i]) 

            if any(isnan.(grz_i)) 
                return true
            end

            push!(grazing, grz_i)
        end
    
        return grazing
    
    else

        #----- for 1N 2P 2Z 2B 2D
        if nz == 2
            K_g = ds["K_g"][2]
            prey = GrM[2,3:end]' .*B[:,1:end]

            for i in range(1, nb)
                gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
                grz_i = (gb_i .* Z[:,2] .* GrM[2,i+2]' ./ prey[:,i])
                push!(grazing, grz_i)
            end

        #----- for 1N 1P 3Z 2B 2D
        elseif nz == 3
            K_g = [ds["K_g"][2], ds["K_g"][3]]
            prey = 1.0 .*B[:,1:2]

            for i in range(1, nb)
                gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g[i])
                grz_i = (gb_i .* Z[:,i+1] .* 1.0 ./ prey[:,i]) 
                push!(grazing, grz_i)
            end
        else
        end
    end

    return grazing

end


function b_loss(mortality, grazing)

    nb = get_size([B])

    if nb == 13
        loss = zeros(Float64, 89, nb)
        for i in range(1, nb)
            loss[:, i] = mortality[i] .+ grazing[:, i]
        end
    else
        loss = Any[]
        for i in range(1, nb)
            push!(loss, mortality[i] .+ grazing[i])
        end
    end
    return loss
end


function RstarB(loss, ds, season)

    season != NaN ? temp_mod = get_temp_mod(season) : temp_mod = ds["temp_fun"][:]

    vmax_ij = ds["vmax_ij"][:]
    Km_ij = ds["Km_ij"][:]
    yield = ds["y_ij"][:]
    temp_mod = get_temp_mod(season)
    II, JJ = get_nonzero_axes(ds["CM"][:])
    RS = Any[]

    nb = get_size([B])
    if nb == 13
        for j = axes(II, 1)
            push!(RS, Km_ij[II[j],JJ[j]] .* loss[:, j] ./ (yield[II[j],JJ[j]] .* vmax_ij[II[j],JJ[j]] .* temp_mod .- loss[:, j]))
        end
    else
        for j = axes(II, 1)
            push!(RS, Km_ij[II[j],JJ[j]] .* loss[j] ./ (yield[II[j],JJ[j]] .* vmax_ij[II[j],JJ[j]] .* temp_mod .- loss[j]))
        end
    end

    for i in range(1, length(RS))
        RS[i] = check_for_negatives(RS[i])
    end

    return RS

end



#-----------------------------------------------------------------------------------
#                                     RSTAR P 
#-----------------------------------------------------------------------------------
function get_rstar_P(P, Z, ds, np, nz, season)
    
    mort_p = p_mortality(P, ds, np)
    grz_p = p_grazing(P, Z, ds, np, nz)
    loss_p = p_loss(mort_p, grz_p, np)
    RstarP = RstarP(loss_p, ds, np, season)

    return RstarP

end


function p_mortality(P, ds, np)

    m_lp = 1e-1
    m_qp = 1e-1
    mort_p = Any[]
    for i in range(1, np)
        push!(mort_p, (m_lp .+ m_qp .* P[:,i]))
    end

    return mort_p

end


function p_grazing(P, Z, ds, np, nz)

    GrM = ds["GrM"][:]
    grazing = Any[]
    g_max = 1.0
    K_g = ds["K_g"]

    if np == 8
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
    
    else

        prey = GrM[1,1:np]' .*P[:,:]
        for i in range(1, np)
            gp_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
            grz_i = (gp_i .* Z[:,1] .* GrM[1,i]') ./ prey[:,i]
            push!(grazing, grz_i)
        end

        return grazing

    end 

end


function p_loss(mortality, grazing, np)

    if np == 8
        loss = zeros(Float64, 89, np)
        for i in range(1, np)
            loss[:, i] = mortality[i] .+ grazing[:, i]
        end
    else
        loss = Any[]
        for i in range(1, np)
            push!(loss, mortality[i] .+ grazing[i])
        end
    end

    return loss

end


function RstarP(loss, ds, np, season=NaN)

    season != NaN ? temp_mod = get_temp_mod(season) : temp_mod = ds["temp_fun"][:]
    umax_ij = ds["umax_ij"][:]
    Kp_ij = ds["Kp_ij"][:]
    temp_mod = get_temp_mod(season)
    RS = Any[]

    if np == 8
        for i in range(1, np)
            push!(RS, Kp_ij[i] .* loss[:, i] ./ (umax_ij[i] .* temp_mod .- loss[:, i]))
        end
    else
        for i in range(1, np)
            push!(RS, (Kp_ij[i] .* loss[i]) ./ (umax_ij[i] .* temp_mod .- loss[i]))
        end
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

# rs = get_rstar_B(bw, zw, winter)