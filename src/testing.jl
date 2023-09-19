using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames, CSV
using SparseArrays, LinearAlgebra


function get_endpoints_ep(vars)

    out = Vector{Any}()

    for v in vars
        append!(out, [v[:,:,end]])
    end

    return out[1], out[2], out[3], out[4], out[5]

end

function get_nonzero_axes(M)

    Cs = sparse(M)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 

function cut_off(ds, n)

    dss = copy(ds)
    co = 10^-6
    for i in range(1, n)
        dss[:, i] .= ifelse.(dss[:, i] .< co, co, dss[:, i])
    end

    return dss

end

function get_zc(ds, H=890)

    dz = 10
    zc = [dz/2:dz:(H-dz/2)]

    return zc
end

function group_BD_competitors(CM)

    Cs = sparse(CM)
    competitors = Any[]

    for (i, row) in enumerate(eachrow(Cs))
        for (j, col) in enumerate(eachcol(Cs))
            if Cs[i, j] > 0
                push!(competitors, [i, j])
            end
        end
    end

    return competitors

end


ds = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230915_22:19_8P6Z13B5D_ep.nc")
II, JJ = get_nonzero_axes(ds["CM"][:])