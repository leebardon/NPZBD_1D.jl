
using DataFrames, NCDatasets, JLD, Printf

function message(v::String, nd::Int64=0, nb::Int64=0, nn::Int64=0, np::Int64=0, nz::Int64=0, fsaven::String="")

    m = Dict(
        "START" => "\n -------------------------------- STARTING PROGRAM ----------------------------------- \n",
        "ST1" => ["Start New Run", "Use Saved Params", "Use Prescribed Model"],
        "ST2" => "\nStart new run, load saved params or run prescribed model?",
        "TM1" => ["2 years (tt=732)", "10 years (tt=3660)", "30 years (tt=10980)", "100 years (tt=36600)"],
        "TM2" => "Select Simulation Runtime:",
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
        "UP2" => "Select max bacteria uptake rate (vmax_i, 1/d):",
        "UPP1" => ["Ordered assignment", "Randomly selected along log range"],
        "UPP2" => "\nSelect max phyto uptake rate (umax_i, 1/d):",
        "MIC" => "\n SETTING MICROBIAL TRAITS \n --------------------------",
        "T1" => ["Tradeoff", "Constant affinity"],
        "TB2" => "Apply bacterial growth rate-affinity tradeoff?",
        "TP2" => "\nApply phyto growth rate-affinity tradeoff?",
        "ENV" => "\n SETTING NUTRIENT SUPPLY \n -------------------------- ",
        "SE2" => "Simulate winter or summer conditions?",
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

    println(message("DN"))
    input = readline()
    nd = parse(Int64, input) 
    nb = microbe_num("BN")
    np = microbe_num("PN")
    nz = microbe_num("ZN")
    nn = 1

    yield = request(message("Y2"), RadioMenu(message("Y1")))
    yield == 1 ? y_i = ones(nd)*0.4 : y_i = rand(nd)*0.5

    supply_weight = request(message("SW2"), RadioMenu(message("SW1")))

    if run_type != 3
        println(message("SUB"))
        uptake = request(message("UP2"), RadioMenu(message("UP1")))
        uptake == 1 ? vmax_i = ordered_uptake_arr(nd) : vmax_i = random_uptake_arr(nd)

        uptake_p = request(message("UPP2"), RadioMenu(message("UPP1")))
        uptake_p == 1 ? umax_i = fill(1., np) : umax_i = random_uptake_arr(np)
    else 
        vmax_i = get_prescribed_params("vmax_i") 
        umax_i = get_prescribed_params("umax_i") 
    end 

    println(message("ENV"))
    season = request(message("SE2"), RadioMenu(message("SE1")))

    @info("User Selections: \n SW = $supply_weight \n B yield = $y_i \n B uptake = $vmax_i \n P uptake = $umax_i \n Season == $season \n")

    return nd, nb, np, nz, nn, y_i, supply_weight, vmax_i, umax_i, season

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


function test_vals(arr)

    e_msg = "\n Nan or inf found in timestep "
    w_msg = "\n Weird values found in timestep "
    check = run_checks(arr)

    return check, arr, e_msg, w_msg

end


function run_checks(vals)

    for x in vals
        if nan_or_inf(x)
            return "e"
        elseif big_or_small(x)
            return "w"
        end
    end

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


function big_or_small(x)

    if typeof(x) == Float64 || typeof(x) == Int64
        if x > 1e10 || x < -1e10 
            return true
        end
    else
        for i in x
            if i > 1e10 || i < -1e10
                return true
            end
        end
    end

    return false

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

    @printf("\n np = %5.0f \n nb = %5.0f \n nz = %5.0f \n nn = %5.0f \n nd = %5.0f \n days = %5.0f \n\n", prms.np, prms.nb, prms.nz, prms.nn, prms.nd, prms.tt)
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

    @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t_id, prms.tt, t_id/prms.tt*100, now())

    return track_n, track_p, track_z, track_b, track_d, track_o, track_time

end


function get_endpoints(vars, ds=nothing)

    endpoints = Vector{Any}()

    if ds !== nothing
        for v in vars
            append!(endpoints, [ds["$v"][:,:,end]])
        end
    else
        for v in vars
            append!(endpoints, [v[:,:,end]])
        end
    end

    return endpoints
end

function set_extinct_to_zero(ds)

    dss = copy(ds)
    ex = 10^-6
    dss .= ifelse.(dss .<= ex, 0.0, dss)

    return dss

end

#------------------------------------------------------------------------
#                             PLOT UTILS 
#------------------------------------------------------------------------

function get_plot_vars()

    bcols = ["cyan3", "darkorange", "indigo", "coral4", "lightcyan4", "magenta2", "thistle", "seagreen4",
    "darkkhaki", "purple", "crimson",  "yellow3", "navajowhite4"]
    dcols = ["blue3", "black", "maroon", "coral", "orange3"]
    pcols = ["olivedrab3", "darkgreen","red4", "cyan4", "gold3", "black", "hotpink2", "wheat2" ]
    ncols = ["blue2"]
    zcols = ["black", "slategray4", "deeppink3", "sienna", "mediumpurple3", "darkseagreen"]
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


# BELOW WAS INSERTED INTO FUNCTIONS TO TRACE NAN AND INF WEIRDNESS - CAUSED BY USING UNDEF TO 
# INITIALISE EMPTY ARRS TO BE LATER USED DURING INTEGRATION


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