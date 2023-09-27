using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames, CSV
using SparseArrays, LinearAlgebra


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