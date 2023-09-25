using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames
using SparseArrays, LinearAlgebra

winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230827_1345.nc")
summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230827_1710.nc")


function rstar_analysis()


end 


function get_rstar_B(B, Z, ds)
    
    mort_b = b_mortality(B, ds)
    grz_b = b_grazing(B, Z, ds)
    loss_b = b_loss(mort_b, grz_b)
    RstarB_ij = star(loss_b, ds)

    return RstarB_ij

end

function b_mortality(B, ds)

    mort_b = Any[]
    for i in range(1, 7)
        push!(mort_b, (ds["m_lb"][i] .+ ds["m_qb"][i] .* B[:,i]) .* B[:,i])
    end

    return mort_b

end

function b_grazing(B, Z, ds)

    GrM = ds["GrM"][:]
    grazing = Any[]
    g_max = 1.0
    K_g = ds["K_g"][3]
    K_g_POM = ds["K_g"][2]

    prey_POM = GrM[2,5]' .*B[:,1]
    gb_POM = g_max .* prey_POM ./ (prey_POM .+ K_g_POM)
    grz_POM = gb_POM .* Z[:,2] .* GrM[2,5]'  .* B[:,1] ./ prey_POM
    push!(grazing, grz_POM)

    prey = GrM[3,6:end]' .*B[:,2:end]
    for i in range(1, 6)
        gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
        grz_i = gb_i .* Z[:,3] .* GrM[3,i+5]'  .* B[:,i+1] ./ prey[:,i]
        push!(grazing, grz_i)
    end

    return grazing

end

function b_loss(mortality, grazing)

    loss = Any[]
    for i in range(1, 7)
        push!(loss, mortality[i] .+ grazing[i])
    end

    return loss
end

function Rstar(loss, ds)

    II, JJ = get_nonzero_axes(ds["CM"][:])
    vmax_ij = ds["vmax_ij"][:]
    Km_ij = ds["Km_ij"][:]

    RS = Any[]
    for j = axes(II, 1)
        push!(RS, Km_ij[II[j],JJ[j]] .* loss[j] ./ (vmax_ij[II[j],JJ[j]] .- loss[j]))
    end

    return RS
end

function get_endpoints(ds, vars)

    out = Vector{Any}()

    for v in vars
        append!(out, [ds["$v"][:,:,end]])
    end

    return out[1], out[2], out[3], out[4], out[5]

end

function get_nonzero_axes(M)

    Cs = sparse(M)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 

nw, pw, zw, bw, dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
ns, ps, zs, bs, ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

rs = get_rstar_B(bw, zw, winter)