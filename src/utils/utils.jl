
using DataFrames, NCDatasets, Printf

function message(v::String, nd::Int64=0, nb::Int64=0, nn::Int64=0, np::Int64=0, nz::Int64=0, fsaven::String="")

    m = Dict(
        "START" => "\n -------------------------------- STARTING PROGRAM ----------------------------------- \n",
        "ST1" => ["Start New Run", "Start Prescribed Model Run", "Continue Run"],
        "ST2" => "\nChoose run type: ",
        "TM1" => ["2 years (days=732)", "10 years (days=3660)", "30 years (days=10980)", "50 years (days=18300)", "100 years (days=36600)"],
        "TM2" => "\nSelect Simulation Runtime:",
        "DR1" => ["Yes", "No"],
        "DR2" => "\n Are you using params from prescribed_darwin.jl?",
        "P1" => ["None (steady state)", "Single pulse at 10 (winter) or 30 (summer) day intervals"],
        "P2" => "Select nutrient pulsing regime: ",
        "DN" => "\nEnter number of detritus pools (nd): ",
        "BN" => "Enter number of bacteria populations (nb): ",
        "PN" => "Enter number of phyto populations (np): ",
        "ZN" => "Enter number of zooplank populations (nz): ",
        "Y1" => ["Equal for all bacteria", "Randomised"],
        "Y2" => "\n Select yield rates for bacteria populations (y_i):",
        "SW1" => ["Prescribed distribution", "Equal distribution", "Lognormal distribution"],
        "SW2" => "\n Select supply weight for OM pools (SW):",
        "SUB" => "\n SETTING SUBSTRATE TRAITS \n -------------------------- ",
        "UP1" => ["Ordered assignment", "Randomly selected along log range"],
        "UP2" => "Select max bacteria uptake rate (umax_i, 1/d):",
        "UPP1" => ["Ordered assignment", "Randomly selected along log range"],
        "UPP2" => "\nSelect max phyto uptake rate (vmax_i, 1/d):",
        "MIC" => "\n SETTING MICROBIAL TRAITS \n --------------------------",
        "T1" => ["Tradeoff", "Constant affinity"],
        "TB2" => "Apply bacterial growth rate-affinity tradeoff?",
        "TP2" => "\nApply phyto growth rate-affinity tradeoff?",
        "ENV" => "\n SETTING NUTRIENT SUPPLY \n -------------------------- ",
        "SE2" => "\nSimulate winter or summer conditions?",
        "SE1" => ["Winter", "Summer"],
        "LVP" => "\nSelect params to load: ",
        "SVP1" => ["Yes", "No"],
        "SVP2" => "\n Save params from this run?",
        "SV" => "Saving to: $fsaven",

    )

    return m["$v"]

end


function microbe_num(MSG)
    println(message(MSG))
    input = readline()
    n = parse(Int64, input)
    return n
end


function user_select(run_type=0)

    if run_type == 1
        println(message("DN"))
        input = readline()
        nd = parse(Int64, input) 
        nb = microbe_num("BN")
        np = microbe_num("PN")
        nz = microbe_num("ZN")
        nn = 1
        println(message("SUB"))
        uptake = request(message("UP2"), RadioMenu(message("UP1")))
        uptake == 1 ? umax_i = ordered_uptake_arr(nd) : umax_i = random_uptake_arr(nd)
        uptake_p = request(message("UPP2"), RadioMenu(message("UPP1")))
        uptake_p == 1 ? vmax_i = fill(1., np) : vmax_i = random_uptake_arr(np)

    elseif run_type == 2
        umax_i = get_prescribed_params("umax_i") 
        vmax_i = get_prescribed_params("vmax_i") 
        nd = length(umax_i)
        np = length(vmax_i)
        nb = length(get_prescribed_params("Fg_b"))
        nz = size(get_prescribed_params("GrM"), 1)
        nn = 1

    end 

    yield = request(message("Y2"), RadioMenu(message("Y1")))
    yield == 1 ? y_i = ones(nd)*0.3 : y_i = rand(nd)*0.5
    supply_weight = request(message("SW2"), RadioMenu(message("SW1")))

    println(message("ENV"))
    pulse = request(message("P2"), RadioMenu(message("P1")))
    season = request(message("SE2"), RadioMenu(message("SE1")))

    @info("User Selections: \n pulse type = $pulse, SW = $supply_weight \n B yield = $y_i \n B uptake = $umax_i \n P uptake = $vmax_i \n Season == $season \n")

    return nd, nb, np, nz, nn, y_i, supply_weight, umax_i, vmax_i, season, pulse

end


function get_previous_params()

    files = readdir("results/outfiles/endpoints")
    f = request("\nSelect output file:", RadioMenu(files))
    ds = NCDataset("results/outfiles/endpoints/$(files[f])")

    H = ds["H"][:][1]
    dz = ds["dz"][:][1]
    nIC, pIC, zIC, bIC, dIC, oIC0 = get_endpoints(["n", "p", "z", "b", "d", "o"], ds)
    oIC = oIC0[:,:]
    nn, np, nz, nb, nd = get_size([nIC, pIC, zIC, bIC, dIC])
    vmax_i = ds["vmax_i"][:]
    vmax_ij = ds["vmax_ij"][:,:]
    Kp_i = ds["Kp_i"][:]
    Kp_ij = ds["Kp_ij"][:,:]
    m_lp = ds["m_lp"][:]
    m_qp = ds["m_qp"][:]
    light = ds["light"][:]
    temp_fun = ds["temp_fun"][:]
    K_I = ds["K_I"][:][1]
    CMp = ds["CMp"][:,:]
    Fg_p = ds["Fg_p"][:]
    umax_i = ds["umax_i"][:]
    umax_ij = ds["umax_ij"][:,:]
    Km_i = ds["Km_i"][:]
    Km_ij = ds["Km_ij"][:,:]
    m_lb = ds["m_lb"][:]
    m_qb = ds["m_qb"][:]
    y_ij = ds["y_ij"][:,:]
    prob_generate_d = ds["prob_generate_d"][:]
    CM = ds["CM"][:,:]
    Fg_b = ds["Fg_b"][:]
    g_max = ds["g_max"][:]
    K_g = ds["K_g"][:]
    γ = ds["γ"][:]
    m_lz = ds["m_lz"][:]
    m_qz = ds["m_qz"][:]
    GrM = ds["GrM"][:,:]
    pen = ds["pen"][:]
    kappa_z = ds["kappa_z"][:]
    wd = ds["wd"][:,:]
    ngrid = ds["ngrid"][:][1]
    pulse = ds["pulse"][:][1]
    e_o = ds["e_o"][:][1]
    yo_ij = ds["yo_ij"][:,:]
    koverh = ds["koverh"][:][1]
    o2_sat = ds["o2_sat"][:][1]
    ml_boxes = ds["ml_boxes"][:][1]
    t_o2relax = ds["t_o2relax"][:][1]
    o2_deep = ds["o2_deep"][:][1]
    season = ds.attrib["Season"]

    return H, dz, nIC, pIC, zIC, bIC, dIC, oIC, nn, np, nz, nb, nd, 
           vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p, 
           umax_i, umax_ij, Km_i, Km_ij, m_lb, m_qb, y_ij, prob_generate_d, CM, Fg_b,
           g_max, K_g, γ, m_lz, m_qz, GrM, pen, kappa_z, wd, ngrid, pulse, e_o, yo_ij,
           koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, season

end


function load_matrix(mtype, nd, nb, nn=0, np=0, nz=0)
    
    if mtype == "CM"
        M = jldopen("results/saved_matrices/$(mtype)_$(nd)d$(nb)b.jdl", "r") do file
            read(file, "A")
        end
    elseif mtype == "CMp"
        M = jldopen("results/saved_matrices/$(mtype)_$(nn)n$(np)p.jdl", "r") do file
            read(file, "A")
        end
    else 
        M = jldopen("results/saved_matrices/$(mtype)_$(np)p$(nb)b$(nz)z.jdl", "r") do file
            read(file, "A")
        end
    end

    return M

end


function save_matrices(M1, M2, M3, nd, nb, nn, np, nz)

    jldopen("results/saved_matrices/CM_$(nd)d$(nb)b.jdl", "w") do file
        write(file, "A", M1)  
    end
    jldopen("results/saved_matrices/CMp_$(nn)n$(np)p.jdl", "w") do file
        write(file, "A", M2)  
    end
    jldopen("results/saved_matrices/GrM_$(np)p$(nb)b$(nz)z.jdl", "w") do file
            write(file, "A", M3)  
    end

end


function get_matrix(Mtype, nd, nb, nn, np, nz) 
    
    if Mtype == "CM"
        if isfile("NPZBD_1D/results/saved_matrices/CM_$(nd)d$(nb)b.jdl")
            M = load_matrix("CM", nd, nb)
        else
            M = build_consumption_matrix(nd, nb)
        end
    elseif Mtype == "CMp"
        if isfile("NPZBD_1D/results/saved_matrices/CMp_$(np)p$(nd)d$(nb)b$(nz)z.jdl")
            M = load_matrix("CMp", nd, nb, nn, np)
        else
            M = build_consumption_matrix(nn, np)
        end  
    else
        if isfile("NPZBD_1D/results/saved_matrices/GrM_$(np)p$(nb)b$(nz)z.jdl")
            M = load_matrix(Mtype, nd, nb, nn, np, nz)
        else
            M = build_grazing_matrix(np, nb, nz)
        end
    end

    return M
end


function check_for_empty_cols(M, n)

    empty_cols = findall(x -> x == 0, sum(M, dims=1))
    s = size(empty_cols)

    if s[1] > 0
        for i in 1:s[1]
            new_col = sprand(Bool, n, 1, 0.5)
            M[:, empty_cols[i][2]] = new_col
        end
        M = check_for_empty_cols(M, n)
    end

    return M

end


function print_info(prms)

    @printf("\n np = %5.0f \n nb = %5.0f \n nz = %5.0f \n nn = %5.0f \n nd = %5.0f \n days = %5.0f \n\n", prms.np, prms.nb, prms.nz, prms.nn, prms.nd, prms.days)
    println("File will be saved as: ", prms.fsaven)
    println("nt = ", prms.nt)

end


function set_logger(launch_time)

    loginfo = string(Dates.format(launch_time, "yyyymmdd_HHMM"), ".log")
    logger = activate_logger(loginfo)

    return logger

end


function activate_logger(loginfo)

    logger = TeeLogger(
        MinLevelLogger(FileLogger("logs/$loginfo"), Logging.Info),
        MinLevelLogger(FileLogger("logs/error.log"), Logging.Warn),
    ); 
    
    global_logger(logger)

    return logger 

end


function log_params(prms)

    @info(
    """Model Params: 
    days:           $(prms.days)
    ts/day (dt):    $(prms.dt)
    total ts (nt):  $(prms.nt)                    
    recorded ts:    $(prms.nt) 
    nn:             $(prms.nn)                    
    np:             $(prms.np)                     
    nz:             $(prms.nz)        
    nb:             $(prms.nb)         
    nd:             $(prms.nd)    
    nIC:            $(prms.nIC[1])                        
    pIC:            $(prms.pIC[1])
    zIC:            $(prms.zIC[1])
    bIC:            $(prms.bIC[1])
    dIC:            $(prms.dIC[1])
    oIC:            $(prms.oIC[1])
    vmax_i:         $(prms.vmax_i)
    vmax_ij:        $(prms.vmax_ij)
    umax_i:         $(prms.umax_i)
    umax_ij:        $(prms.umax_ij)
    Kp_i:           $(prms.Kp_i)
    Kp_ij:          $(prms.Kp_ij)
    Km_i:           $(prms.Km_i)
    Km_ij:          $(prms.Km_ij)
    y_ij:           $(prms.y_ij[1])
    K_g:            $(prms.K_g[1])
    g_max:          $(prms.g_max[1])
    γ:              $(prms.γ[1])
    m_lp:           $(prms.m_lp[1])
    m_qp:           $(prms.m_qp[1])
    m_lb:           $(prms.m_lb[1])
    m_qb:           $(prms.m_qb[1])
    m_lz:           $(prms.m_lz[1])
    m_qz:           $(prms.m_qz[1])
    prob_generate_d: $(prms.prob_generate_d)
    wd:             $(prms.wd)
    pulse? (1=no)   $(prms.pulse)   
    Fg_p:           $(prms.Fg_p)
    Fg_b:           $(prms.Fg_b)
    light half sat: $(prms.K_I)
    """) 
    
    print_info(prms)


end


function update_tracking_arrs(track_n, track_p, track_z, track_b, track_d, track_o, track_time, ntemp, ptemp, ztemp, btemp, dtemp, otemp, t, trec, prms)

    j = Int(t÷trec + 1)
    t_id = t.*prms.dt
    track_p[:,:,j] .= ptemp
    track_b[:,:,j] .= btemp 
    track_z[:,:,j] .= ztemp 
    track_n[:,:,j] .= ntemp 
    track_d[:,:,j] .= dtemp
    track_o[:,:,j] .= otemp
    track_time[j] = t_id 

    @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t_id, prms.days, t_id/prms.days*100, now())

    return track_n, track_p, track_z, track_b, track_d, track_o, track_time

end


function nan_or_inf(x)

    if typeof(x) == Float64 || typeof(x) == Int64
        if isnan(x) || !isfinite(x)
            return true
        end
    else
        if any(isnan.(x)) || any(isinf.(x))
            return true
        end
    end

    return false

end


function get_endpoints(vars, ds=nothing)

    endpoints = Vector{Any}()

    if ds !== nothing
        for v in vars
            if v != "o"
                append!(endpoints, [ds["$v"][:,:,end]])
            else
                append!(endpoints, [ds["$v"][:,end]])
            end
        end
    else
        for v in vars
            if v != "o"
                append!(endpoints, [v[:,:,end]])
            else
                append!(endpoints, [v[:,end]])
            end
        end
    end

    return endpoints
end


function get_final_year(ds, vars)
    # where year is 12 * 30 days
    final_yr = Vector{Any}()

    for v in vars
        if v != "o"
            append!(final_yr, [ds[v][:, :, end-7199:end]])
        else
            append!(final_yr, [ds[v][:, end-7199:end]])
        end
    end

    return final_yr

end


function get_zc(H)

    dz = 10
    zc = [dz/2 : dz : H - dz/2;] 

    return zc

end


function set_extinct_to_zero(ds)

    dss = copy(ds)
    ex = 10^-6
    dss .= ifelse.(dss .<= ex, 0.0, dss)

    return dss

end


function mean_over_time(state_vars, ds, season)

    if season == "Winter"
        pulse_freq = 10
        means_over_time = get_cycle_means(state_vars, pulse_freq, ds)
    else
        pulse_freq = 30
        means_over_time = get_cycle_means(state_vars, pulse_freq, ds)
    end

    return means_over_time
    
end


function get_cycle_means(vars, pulse_freq, ds)

    num_days = pulse_freq * 100
    ts = num_days * 20

    cycle_mean = Vector{Any}()

    for v in vars
        if v != "o"
            period = ds["$v"][:,:,end-ts:end]
            period_mean = dropdims(mean(period, dims=3), dims=3)
            append!(cycle_mean, [period_mean])
        else
            append!(cycle_mean, [mean(ds["$v"][:,end-ts:end], dims=2)])
        end
    end


    return cycle_mean

end




function check_for_negatives(RS)

    for i in eachindex(RS)
        for j in eachindex(RS[i])
            RS[i][j] = ifelse(RS[i][j] < 0 || RS[i][j] > 15, NaN, RS[i][j])
        end
    end

    return RS

end


function get_temp_mod(ds)
    
    try
        temp_mod = ds["temp_fun"][:]
        return temp_mod
    catch
        rgx = r"(?<=/)([A-Za-z][A-Za-z])(?=\d)"
        season = match(rgx, fsaven)
        if season[1] == "Wi"
            temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/win_temp_mod.csv", DataFrame)
        else
            temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/sum_temp_mod.csv", DataFrame)
        end
        return Matrix(temp_mod)
    end

end


function get_light()
    return CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/light/light.csv", DataFrame, header=false)
end


#------------------------------------------------------------------------
#                             PLOT UTILS 
#------------------------------------------------------------------------

function get_plot_vars()

    bcols = ["cyan3", "darkorange", "indigo", "coral4", "lightcyan4", "magenta2", "thistle", "seagreen4",
            "darkkhaki", "purple", "black",  "yellow3", "navajowhite4",  "coral4", "orange2", "orangered4", "yellow3", 
            "lightyellow4", "goldenrod4", "slateblue4", "mediumpurple3"]
    dcols = ["blue3", "black", "maroon", "coral", "orange3", "silver", "magenta3", "cyan3", "seagreen", "mediumpurple3"]
    pcols = ["hotpink2", "darkgreen","red4", "cyan4", "gold3", "black", "brown", "wheat2", "mediumpurple3", "darkseagreen" ]
    ncols = ["blue2"]
    zcols = ["black", "slategray4", "deeppink3", "sienna", "mediumpurple3", "darkseagreen", "snow4", 
            "silver", "salmon", "coral4", "orange2", "orangered4", "yellow3", "lightyellow4", "goldenrod4",
            "chartreuse", "lightseagreen", "blueviolet", "slateblue4", "magenta4" ]
    ab = 0.8
    ab_ext = 0.8
    ls = 4
    lfs = 9
    lg = :bottomright
    
    return bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg

end


function get_size(arr)

    out = Vector{Int}()
    
    for a in arr
        append!(out, size(a, 2))
    end

    return out

end


function get_nonzero_axes(Mat)

    Cs = sparse(Mat)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 


function group_interactions(Cs, n)

    interactions = Any[]
    for (i, row) in enumerate(eachrow(Cs))
        for (j, col) in enumerate(eachcol(Cs))
            if Cs[i, j] > 0
                push!(interactions, [i, j])
            end
        end
    end

    return get_interaction_dict(interactions, n) 

end


function get_interaction_dict(interactions, n)

    out = Dict(name => Any[] for name in collect(1:1:n))
    for i in interactions
        for j in keys(out) 
            if i[1] == j 
                push!(out[j], i[2]) 
            end 
        end 
    end

    return out

end


# BELOW WAS INSERTED INTO FUNCTIONS TO TRACE NAN AND INF WEIRDNESS - CAUSED BY USING UNDEF TO 
# INITIALISE EMPTY ARRS TO BE LATER USED DURING INTEGRATION

# function test_vals(arr)

#     e_msg = "\n Nan or inf found in timestep "
#     w_msg = "\n Weird values found in timestep "
#     check = run_checks(arr)

#     return check, arr, e_msg, w_msg

# end


# function run_checks(vals)

#     for x in vals
#         if nan_or_inf(x)
#             return "e"
#         elseif big_or_small(x)
#             return "w"
#         end
#     end

# end


# function big_or_small(x)

#     if typeof(x) == Float64 || typeof(x) == Int64
#         if x > 1e10 || x < -1e10 
#             return true
#         end
#     else
#         for i in x
#             if i > 1e10 || i < -1e10
#                 return true
#             end
#         end
#     end

#     return false

# end

# check, data, e_msg, w_msg = test_vals([dNdt, dPdt, dZdt, dBdt, dDdt])
# if check=="e"
#     @error("$e_msg $t at j=$j: \n $data \n")
# elseif check=="w"
#     print("warn")
#     @error("$w_msg $t at j=$j: \n $data \n")
# end

#         #TODO find a better way to do this so can traceback source of fault and not need to repeat code blocks
#         #NOTE this is probably something that can be improved with proper unit testing
#         check, data, e_msg, w_msg = test_vals([uptake, mort, dNdt, dPdt, d_gain_total])
#         if check=="e"
#             @error("$e_msg $t at i=$i: \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t at i=$i: \n $data \n ")
#         end


#         check, data, e_msg, w_msg = test_vals([uptake, dDdt, dBdt, dNdt])
#         if check=="e"
#             @error("$e_msg $t at j=$j: \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t at j=$j: \n $data \n")
#         end


#         check, data, e_msg, w_msg = test_vals([prey, g, dZdt, dNdt, dPdt])
#         if check=="e"
#             @error("$e_msg $t at k=$k: \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t at k=$k: \n $data \n")
#         end  


#         check, data, e_msg, w_msg = test_vals([prey, g, dZdt, dNdt, dBdt])
#         if check=="e"
#             @error("$e_msg $t at k=$k: \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t at k=$k: \n $data \n")
#         end


#         check, data, e_msg, w_msg = test_vals([bmort, dBdt, d_gain_total])
#         if check=="e"
#             @error("$e_msg $t : \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t : \n $data \n")
#         end


#         check, data, e_msg, w_msg = test_vals([zmort, dZdt, d_gain_total])
#         if check=="e"
#             @error("$e_msg $t : \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t : \n $data \n")
#         end


#         check, data, e_msg, w_msg = test_vals([dDdt])
#         if check=="e"
#             @error("$e_msg $t : \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t : \n $data \n")
#         end