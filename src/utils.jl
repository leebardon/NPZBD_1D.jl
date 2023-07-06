
function message(v::String, nd::Int64=0, nb::Int64=0, nn::Int64=0, np::Int64=0, nz::Int64=0, fsaven::String="")

    m = Dict(
        "START" => "\n -------------------------------- STARTING PROGRAM ----------------------------------- \n",
        "ST1" => ["Start New Run", "Use Saved Params"],
        "ST2" => "\nStart new run, or load saved params?",
        "TM1" => ["1 year (tt=366)", "10 years (tt=3660)", "30 years (tt=10980)", "100 years (tt=36600)"],
        "TM2" => "Select Simulation Runtime:",
        "DN" => "\nEnter number of detritus pools (nd): ",
        "BN" => "Enter number of bacteria populations (nb): ",
        "PN" => "Enter number of phyto populations (np): ",
        "ZN" => "Enter number of zooplank populations (nz): ",
        "Y1" => ["Equal for all bacteria", "Randomised"],
        "Y2" => "\n Select yield rates for bacteria populations (y_i):",
        "SW1" => ["Equal distribution", "Lognormal distribution"],
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


function user_select()

    println(message("DN"))
    input = readline()
    nd = parse(Int64, input) 

    nb = microbe_num("BN")
    if !iseven(nb)
        println("\n !!! PLEASE ENTER EVEN NUMBER !!! \n\n")
        nb = microbe_num("BN")
    end

    np = microbe_num("PN")
    if !iseven(np)
        println("\n !!! PLEASE ENTER EVEN NUMBER !!! \n\n")
        np = microbe_num("PN")
    end

    nz = microbe_num("ZN")
    nn = 1

    yield = request(message("Y2"), RadioMenu(message("Y1")))
    supply_weight = request(message("SW2"), RadioMenu(message("SW1")))
    println(message("SUB"))
    uptake = request(message("UP2"), RadioMenu(message("UP1")))
    uptake_p = request(message("UPP2"), RadioMenu(message("UPP1")))
    println(message("ENV"))
    season = request(message("SE2"), RadioMenu(message("SE1")))

    @info("User Selections: \n SW = $supply_weight \n B yield = $yield \n B uptake = $uptake \n P uptake = $uptake_p \n Season == $season \n")

    return nd, nb, np, nz, nn, yield, supply_weight, uptake, uptake_p, season

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


function print_info(prms)

    @printf("\n np = %5.0f \n nb = %5.0f \n nz = %5.0f \n nn = %5.0f \n nd = %5.0f \n days = %5.0f \n\n", prms.np, prms.nb, prms.nz, prms.nn, prms.nd, prms.tt)
    println("File will be saved as: ", prms.fsaven)
    println("nt = ", prms.nt)

end


function set_savefiles(launch_time, years)

    fsave = "results/outfiles/out_$(years)y"
    fsaven = string(fsave, "_", Dates.format(launch_time, "yyyymmdd_HHMM"), ".nc")

    loginfo = replace(fsaven, "out_$(years)y_" => "", ".nc" => ".log", "results/outfiles/" => "")
    logger = activate_logger(loginfo)

    return fsaven, logger

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


function savetoNC(p, b, z, n, d, o, timet, v, uptake, tst, tfn, prms, season_num)

    outdir = "/home/lee/Dropbox/Development/NPZBD_1D/"
    path = joinpath(outdir, prms.fsaven) 
    println("\nSaving output file to: ", path)

    season_num == 1 ? season = "winter" : season = "summer"

    f = NCDataset(path, "c") 

    # define the dim of p, b, z, n, d
    defDim(f,"np", prms.np)
    defDim(f,"nb",prms.nb)
    defDim(f,"nz",prms.nz)
    defDim(f,"nn",prms.nn)
    defDim(f,"nd",prms.nd)

    # define the dim of the depth
    defDim(f,"ndepth",prms.ngrid)
    defDim(f,"ndepth1",prms.ngrid+1)

    # define the dim of the time length
    nrec1 = Int(prms.nrec+1) #bc i added time 0
    nprey = prms.np + prms.nb
    
    defDim(f,"nrec",nrec1)
    defDim(f,"nprey",nprey)
   
    # info
    f.attrib["title"] = "NPZBD 0D model i/o"
    f.attrib["Start time"] = string(tst)
    f.attrib["End time"] = string(tfn)
    f.attrib["Run time"] = string(tfn - tst) 
    f.attrib["Season"] = season

    # simulated results
    w = defVar(f,"p",Float64,("ndepth" ,"np","nrec"))
    w[:,:,:] = p
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"b",Float64,("ndepth" ,"nb","nrec"))
    w[:,:,:] = b
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"z",Float64,("ndepth" ,"nz","nrec"))
    w[:,:,:] = z
    w.attrib["units"] = "mmol/m3 C biomass"
    
    w = defVar(f,"n",Float64,("ndepth" ,"nn","nrec"))
    w[:,:,:] = n
    w.attrib["units"] = "mmol/m3 C OM"

    w = defVar(f,"d",Float64,("ndepth" ,"nd","nrec"))
    w[:,:,:] = d
    w.attrib["units"] = "mmol/m3 C OM"
    
    w = defVar(f,"o",Float64,("ndepth" ,"nrec"))
    w[:,:] = o
    w.attrib["units"] = "mmol/m3 O2"

    # --------------------------------------------------
    
    w = defVar(f,"pIC",Float64,("ndepth","np"))
    w[:,:] = prms.pIC
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"bIC",Float64,("ndepth","nb"))
    w[:,:] = prms.bIC
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"zIC",Float64,("ndepth","nz"))
    w[:,:] = prms.zIC
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"nIC",Float64,("ndepth","nn"))
    w[:,:] = prms.nIC
    w.attrib["units"] = "mmol/m3 C OM"

    w = defVar(f,"dIC",Float64,("ndepth","nd"))
    w[:,:] = prms.dIC
    w.attrib["units"] = "mmol/m3 C OM"

    # --------------------------------------------------
    
    w = defVar(f, "timet", Float64, ("nrec",))
    w[:] = timet
    w.attrib["units"] = "days"

    w = defVar(f, "H", Int, ())
    w[:] = prms.H
    w.attrib["units"] = "m; total height"

    w = defVar(f, "dz", Int, ())
    w[:] = prms.dz
    w.attrib["units"] = "m; box height"

    w = defVar(f, "kappa_z", Float64, ("ndepth1",))
    w[:] = prms.kappa_z
    w.attrib["units"] = "vertical water velocity"
    
    w = defVar(f, "wd", Float64, ("ndepth1","nd"))
    w[:] = prms.wd
    w.attrib["units"] = "sinking rate"

    # --------------------------------------------------

    w = defVar(f,"uptake",Float64,("nd","nb"))
    w[:,:,:] = uptake
    w.attrib["units"] = "mmol/m3 C per d; uptake matrix"

    w = defVar(f,"SW",Float64,("nd",))
    w[:,:] = prms.prob_generate_d 
    w.attrib["units"] = "Ind C supply weight: probability"
    
    # w = defVar(f,"SW_all",Float64,("nd","nrec"))
    # w[:,:] = SW_all
    # w.attrib["units"] = "Ind C supply weight: over time"

    w = defVar(f,"pen",Float64,("nb",))
    w[:,:] = prms.pen
    w.attrib["units"] = "penalty"
    
    w = defVar(f, "umax_i", Float64, ("np",))
    w[:] = prms.umax_i
    w.attrib["units"] = "m3/mmol/d; max growth rate of p"

    w = defVar(f, "umax_ij", Float64, ("nn", "np"))
    w[:,:] = prms.umax_ij
    w.attrib["units"] = "per n; max uptake rate"

    w = defVar(f, "Kp_ij", Float64, ("nn", "np"))
    w[:] = prms.Kp_ij
    w.attrib["units"] = "mmol/m3; half-sat"

    w = defVar(f,"CM",Float64,("nd","nb"))
    w[:,:] = prms.CM
    w.attrib["units"] = "Consumption Matrix (nb x nd)"

    w = defVar(f,"CMp",Float64,("nn","np"))
    w[:,:] = prms.CMp
    w.attrib["units"] = "Consumption Matrix (np x nn)"

    w = defVar(f,"GrM",Float64,("nz","nprey"))
    w[:,:] = prms.GrM
    w.attrib["units"] = "Grazing Matrix"

    w = defVar(f, "y_ij", Float64, ("nd","nb"))
    w[:,:] = prms.y_ij
    w.attrib["units"] = "per d; max yield rate"

    w = defVar(f, "vmax_i", Float64, ("nd",))
    w[:,:] = prms.vmax_i
    w.attrib["units"] = "per d; max uptake rate"
    
    w = defVar(f, "vmax_ij", Float64, ("nd", "nb"))
    w[:,:] = prms.vmax_ij
    w.attrib["units"] = "per d; max uptake rate"
    
    w = defVar(f, "Km_ij", Float64, ("nd", "nb"))
    w[:] = prms.Km_ij
    w.attrib["units"] = "mmol/m3; half-sat"
    
    w = defVar(f, "m_lb", Float64, ("nb",))
    w[:,:] = prms.m_lb
    w.attrib["units"] = "m3/mmol; linear death rate of b"

    w = defVar(f, "m_qb", Float64, ("nb",))
    w[:,:] = prms.m_qb
    w.attrib["units"] = "m3/mmol; quadratic death rate of b"
    
    w = defVar(f, "K_g", Float64, ("nz",))
    w[:] = prms.K_g
    w.attrib["units"] = "m3/mmol; half-sat rate of z"
    
    w = defVar(f, "γ", Float64, ("nz",))
    w[:] = prms.γ
    w.attrib["units"] = "fraction of digestion of z"
    
    w = defVar(f, "m_lz", Float64, ("nz",))
    w[:,:] = prms.m_lz
    w.attrib["units"] = "m3/mmol; linear death rate of z"

    w = defVar(f, "m_qz", Float64, ("nz",))
    w[:,:] = prms.m_qz
    w.attrib["units"] = "m3/mmol; quadratic death rate of z"
    
    close(f)

end

function define_dims(ds, prms, nrec1)

    vars = Dict("nb" => prms.nb, "nd" => prms.nd, "nrec" => nrec1)

    for (k, v) in x
        defDim(ds, k, v)
    end

    return ds

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