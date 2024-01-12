using NCDatasets, SparseArrays, LinearAlgebra, Statistics
using CSV, DataFrames, LaTeXStrings, Plots

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/plotting/productivity_plots.jl")


function get_productivity(fsaven)

    ds = NCDataset(fsaven)
    ds.attrib["Season"] == "winter" ? season_num = 1 : season_num = 2

    if ds["pulse"][:] == 1
        N, P, Z, B, D = get_endpoints(["n", "p", "z", "b", "d"], ds)
    else
        N, P, Z, B, D = mean_over_time(["n", "p", "z", "b", "d"], ds, season_num)
    end

    BP = bacterial_productivity(B, D, ds)
    mean_Leu, mean_Thy = get_BP_spot(ds.attrib["Season"])

    fig, dir, filename = plot_seasonal_BP(BP, mean_Leu, mean_Thy, ds, fsaven, season_num)
    savefig(fig,"$(dir)/$(filename).png")

end


function bacterial_productivity(B, D, ds)

    B[B .< 1e-6] .= 0
    D[D .< 1e-6] .= 0

    CM = ds["CM"][:,:]
    umax_ij = ds["umax_ij"][:,:]
    Km_ij = ds["Km_ij"][:,:]
    temp = ds["temp_fun"][:]

    II, JJ = get_nonzero_axes(CM)
    BP = zeros(89, 13)

    for j = axes(II, 1)
        mu = temp .* umax_ij[II[j],JJ[j]] .* D[:,II[j]] ./ (D[:,II[j]] .+ Km_ij[II[j],JJ[j]])
        BP[:,JJ[j]] = B[:,JJ[j]] .* mu
    end

    BP_out = convert_units(BP)

    return BP

end



function get_BP_spot(season)

    csv_path = "data/spot_data/seasonal_means/seasonal_means.csv"
    df = DataFrame(CSV.File(csv_path))

    leu_df, thy_df = select_cols(df)
    mean_leu, mean_thy = filter_by_season(leu_df, thy_df, season)

    return mean_leu, mean_thy

end


function select_cols(df)

    leu_df = select(df[:,:], :depth, :mean_Leu, :sd_Leu, :season)
    thy_df = select(df[:,:], :depth, :mean_Thy, :sd_Thy, :season)

    return leu_df, thy_df
end


function filter_by_season(leu_df, thy_df, seson)

    leu = filter(row -> row.season == seson, leu_df)
    thy = filter(row -> row.season == seson, thy_df)

    return leu, thy

end 


function convert_units(BP)

    # Estimate of cell nitrogen quota -> 1 fmol N/cell = 1e^-12 mmol N/cell
    # 1 mmol N/1e^-12 cells = 1e^12 cells per mol N ; 1 m3 = 1e^6 mL
    # 1 mmol N/m3 == 1e^12 / 1e^6 = 1e^6 cells/mL

    return BP .*= 1e6


end

# fsaven = "results/outfiles/Wi50y_240103_14:18_6P3Z13B8D.nc" # win steady
# fsaven = "results/outfiles/Wi50y_240108_20:47_6P3Z13B8D.nc" # win pulsed
# fsaven="results/outfiles/Su50y_240104_14:04_6P3Z13B8D.nc" # sum steady
# fsaven="results/outfiles/Su50y_240109_00:03_6P3Z13B8D.nc" #sum pulsed
# get_productivity(fsaven)


# function get_spot_productivity()

#     csv_path = "data/spot_data/SPOT_CTD_Nutrients_BP_5Depths.csv"
#     df = DataFrame(CSV.File(csv_path))

#     leu_BP, thy_BP = get_BP(df)

#     return leu, thy

# end

#TODO group all measurements by depth bin, get average DCM depth in each season
# then get average BP for each depthbin for each season
# use average DCM value for plotting
#NOTE remembered I had already processed this with earlier python work! Added csv's to spot data folder
# function get_spot_productivity()

#     csv_path = "data/spot_data/SPOT_CTD_Nutrients_BP_5Depths.csv"
#     df = DataFrame(CSV.File(csv_path))

#     leu_BP, thy_BP = get_BP(df)

#     return leu, thy

# end


# function get_BP(df)

#     leu_df, thy_df = select_cols(df)

#     mean_Leu, mean_Thy = replace_missing(leu_df), replace_missing(thy_df)
#     leu_, thy_ = to_float(mean_Leu), to_float(mean_Thy)
#     leu, thy = to_cells_per_m3(leu_), to_cells_per_m3(thy_)

#     return leu, thy
# end


# function replace_missing(df)

#     out = ifelse.(df .== "No data", missing, df)

#     return out
# end


# function select_cols(df)

#     leu_df = select(df[3:end,:], :depth, :mean_Leu, :sd_Leu)
#     thy_df = select(df[3:end,:], :depth, :mean_Thy, :sd_Thy)

#     return leu_df, thy_df
# end


# function to_float(df)

#     return passmissing(parse).(Float64, df)
# end


# function to_cells_per_m3(df)

#     df[:, 2:3] = df[:, 2:3] .* 1e6

#     return df
# end


# function plot_spot_vs_model(spot, model)


# end


# function make_chunks(depth, n)

#     c = length(depth) ÷ n
#     return [depth[1+c*i:(i == n-1 ? end : c*i+c)] for i = 0:n-1]

# end








# fsaven = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_231019_17:22_10P3Z18B8D.nc"
# get_npp(fsaven, 30)
