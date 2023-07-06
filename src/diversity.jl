using CSV, DataFrames
    

function load_spot_data()

    path = "NPZBD_1D/data/spot_diversity"
    diversity = DataFrame(CSV.File("$(path)/spot_shannon_means.csv"))

    return diversity

end


spot_diversity = load_spot_data()