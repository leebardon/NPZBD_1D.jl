using AutomaticDocstrings
import REPL
using REPL.TerminalMenus
using Printf, NCDatasets
using Dates
using SparseArrays



function message(v::String, nd::Int64=0, nb::Int64=0, np::Int64=0, nz::Int64=0, fsaven::String="")

    m = Dict(
        "START" => "\n -------------------------------- STARTING PROGRAM ----------------------------------- \n",
        "DEF" => "Default values: \n tt = 5 \n nrec = 100 \n nn = 1 \n np = 1 \n nz = 2 \n nb = 6 \n nd = 3 \n y_i = conts. 0.3 \n SW = const. 1.0 \n vmax_i = ordered assignment \n ",
        "DF1" => ["Use defaults", "Select custom params"],
        "DF2" => "Proceed with defaults or select custom params?",
        "T" => "\n Enter simulation run time (tt - number of days): ",
        "REC" => "Enter number of timepoints to record (nrec): ",
        "OM" => "Enter number of organic matter pools (nd): ",
        "BA" => "Enter number of bacteria populations (nb): ",
        "NPZ" => """\n 
        -----------------------------------
        Currently running with default of: 
        1 phyto (np) \n  2 zoo (nz) \n  1 inorganic matter pool (nn) 
        -----------------------------------""",
        "CM" => "\n >> Building consumption matrix with $nd OM pools and $nb bac populations << ",
        "GM" => "\n >> Building grazing matrix with $np phy, $nb bac and $nz zoo populations << ",
        "Y1" => ["Equal for all bacteria", "Randomised"],
        "Y2" => "\n Select yield rates for bacteria populations (y_i):",
        "SW1" => ["Equal distribution", "Randomised"],
        "SW2" => "\n Select supply weight for OM pools (SW):",
        "SUB" => "\n SETTING SUBSTRATE TRAITS \n -------------------------- ",
        "UP1" => ["Ordered assignment", "Randomly selected along log range"],
        "UP2" => "Select max uptake rate (vmax_i, 1/d):",
        "MIC" => "\n SETTING MICROBIAL TRAITS \n --------------------------",
        "GA1" => ["Tradeoff", "Constant affinity"],
        "GA2" => "Apply bacterial growth rate-affinity tradeoff?",
        "SV" => "Saving to: $fsaven",

    )

    return m["$v"]

end


function get_defaults()
    tt = 5
    nrec = 100 
    nn = 1
    np = 1
    nz = 2
    nb = 6
    nd = 3
    y_i = ones(nd)*0.3 
    SW = ones(nd) 
    vmax_i = ordered_vmax(nd)

    return tt, nrec, nd, nb, np, nz, nn, y_i, SW, vmax_i

end


function bacteria_num()
    println(message("BA"))
    input = readline()
    nb = parse(Int64, input)
    return nb
end


function user_select()

    println(message("T"))
    input = readline()
    tt = parse(Int64, input) 

    println(message("REC"))
    input = readline()
    nrec = parse(Int64, input) 

    println(message("OM"))
    input = readline()
    nd = parse(Int64, input) 

    nb = bacteria_num()
    if !iseven(nb)
        println("\n !!! PLEASE SELECT EVEN NUMBER OF BACTERIA !!! \n\n")
        nb = bacteria_num()
    end

    println(message("NPZ"))
    np = 1
    nz = 2
    nn = 1

    yield = request(message("Y2"), RadioMenu(message("Y1")))
    OM_supply_weight = request(message("SW2"), RadioMenu(message("SW1")))
    println(message("SUB"))
    uptake = request(message("UP2"), RadioMenu(message("UP1")))
    input = readline()

    return tt, nrec, nd, nb, np, nz, nn, yield, OM_supply_weight, uptake

end


function load_matrix(mtype, nd, nb, np=0, nz=0)
    
    if nz == 0       
        M = jldopen("results/matrices/$(mtype)_$(nd)dx$(nb)b.jdl", "r") do file
            read(file, "A")
        end
    else
        M = jldopen("results/matrices/$(mtype)_($(np)p+$(nb)b)x$(nz)z.jdl", "r") do file
            read(file, "A")
        end
    end

    return M

end


function save_matrix(M, mtype, nd, nb, np=0, nz=0)

    if nz == 0
        jldopen("results/matrices/$(mtype)_$(nd)dx$(nb)b.jdl", "w") do file
            write(file, "A", M)  
        end
    else
        jldopen("results/matrices/$(mtype)_($(np)p+$(nb)b)x$(nz)z.jdl", "w") do file
            write(file, "A", M)  
        end
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


function print_info(start_time, prms, nt)

    @printf("\n np = %5.0f \n nb = %5.0f \n nz = %5.0f \n nn = %5.0f \n nd = %5.0f \n T = %5.0f \n\n", params.np, params.nb, params.nz, params.nn,  params.nd, params.tt)
    open("jlog.txt","w") do f
        write(f,@sprintf("np = %5.0f, nb = %5.0f, nz = %5.0f, nn = %5.0f, nd = %5.0f, T = %5.0f \n", params.np, params.nb, params.nz, params.nn,  params.nd, params.tt))
    end

    fsaven = string(prms.fsave,"_", Dates.format(start_time, "yyyymmdd"), ".nc")
    if isfile(fsaven)
        fsaven = string(prms.fsave, "_", Dates.format(start_time, "yyyymmdd_HHMM"), ".nc")
    end

    println("Starting time: $start_time, \nFile will be saved as: $fsaven")

    println("\nConsumption Matrix (CM):")
    display("text/plain", CM)

    println("\nGrazing Matrix (GrM):")
    display("text/plain", GrM)

    println("nt = ", nt)

    return fsaven

end

function update_tracking_arrs(track_n, track_p, track_z, track_b, track_d, track_o, track_time, ntemp, ptemp, ztemp, btemp, dtemp, otemp, t, trec, prms)

    j = Int(t÷trec + 1)
    t_id = t.*dt
    track_p[:,:,j] .= ptemp
    track_b[:,:,j] .= btemp 
    track_z[:,:,j] .= ztemp 
    track_n[:,:,j] .= ntemp 
    track_d[:,:,j] .= dtemp
    track_o[:,:,j] .= otemp
    track_time[j] = t_id 


    @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t_id, prms.tt, t_id/prms.tt*100, now())
    open("jlog.txt","a") do f
        write(f,@sprintf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t_id, prms.tt, t_id/prms.tt*100, now()))
    end

    return track_n, track_p, track_z, track_b, track_d, track_o, track_time

end


function savetoNC(fsaven, p, b, z, n, d, o, timet, v, uptake, tst, tfn, params)

    println("Saving to: ",fsaven)

    f = NCDataset(fsaven, "c") #c for create

    # define the dim of p, b, z, n, d
    defDim(f,"np", params.np)
    defDim(f,"nb",params.nb)
    defDim(f,"nz",params.nz)
    defDim(f,"nn",params.nn)
    defDim(f,"nd",params.nd)

    # define the dim of the depth
    defDim(f,"ndepth",params.ngrid)
    defDim(f,"ndepth1",params.ngrid+1)

    # define the dim of the time length
    nrec1 = Int(params.nrec+1) #bc i added time 0
    nprey = params.np + params.nb
    
    defDim(f,"nrec",nrec1)
    defDim(f,"nprey",nprey)
   
    # info
    f.attrib["title"] = "NPZBD 0D model i/o"
    f.attrib["Start time"] = string(tst)
    f.attrib["End time"] = string(tfn)
    f.attrib["Run time"] = string(tfn - tst) 

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
    w[:,:] = params.pIC
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"bIC",Float64,("ndepth","nb"))
    w[:,:] = params.bIC
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"zIC",Float64,("ndepth","nz"))
    w[:,:] = params.zIC
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"nIC",Float64,("ndepth","nn"))
    w[:,:] = params.nIC
    w.attrib["units"] = "mmol/m3 C OM"

    w = defVar(f,"dIC",Float64,("ndepth","nd"))
    w[:,:] = params.dIC
    w.attrib["units"] = "mmol/m3 C OM"

    # --------------------------------------------------
    
    w = defVar(f, "timet", Float64, ("nrec",))
    w[:] = timet
    w.attrib["units"] = "days"

    w = defVar(f, "H", Int, ())
    w[:] = params.H
    w.attrib["units"] = "m; total height"

    w = defVar(f, "dz", Int, ())
    w[:] = params.dz
    w.attrib["units"] = "m; box height"

    w = defVar(f, "kappa_z", Float64, ("ndepth1",))
    w[:] = params.kappa_z
    w.attrib["units"] = "vertical water velocity"
    
    w = defVar(f, "wd", Float64, ("ndepth1","nd"))
    w[:] = params.wd
    w.attrib["units"] = "sinking rate"

    # --------------------------------------------------

    w = defVar(f,"uptake",Float64,("nd","nb"))
    w[:,:,:] = uptake
    w.attrib["units"] = "mmol/m3 C per d; uptake matrix"

    w = defVar(f,"SW",Float64,("nd",))
    w[:,:] = params.SW 
    w.attrib["units"] = "Ind C supply weight: probability"
    
    # w = defVar(f,"SW_all",Float64,("nd","nrec"))
    # w[:,:] = SW_all
    # w.attrib["units"] = "Ind C supply weight: over time"

    w = defVar(f,"pen",Float64,("nb",))
    w[:,:] = params.pen
    w.attrib["units"] = "penalty"
    
    w = defVar(f, "umax_p", Float64, ("np",))
    w[:] = params.umax_p
    w.attrib["units"] = "m3/mmol/d; max growth rate of p"

    w = defVar(f, "K_n", Float64, ("np",))
    w[:] = params.K_n
    w.attrib["units"] = "m3/mmol; half-sat rate of p"

    w = defVar(f,"CM",Float64,("nd","nb"))
    w[:,:] = params.CM
    w.attrib["units"] = "Consumption Matrix"

    w = defVar(f,"GrM",Float64,("nz","nprey"))
    w[:,:] = params.GrM
    w.attrib["units"] = "Grazing Matrix"

    w = defVar(f, "y_ij", Float64, ("nd","nb"))
    w[:,:] = params.y_ij
    w.attrib["units"] = "per d; max yield rate"

    w = defVar(f, "vmax_i", Float64, ("nd",))
    w[:,:] = params.vmax_i
    w.attrib["units"] = "per d; max uptake rate"
    
    w = defVar(f, "vmax_ij", Float64, ("nd", "nb"))
    w[:,:] = params.vmax_ij
    w.attrib["units"] = "per d; max uptake rate"
    
    w = defVar(f, "Km_ij", Float64, ("nd", "nb"))
    w[:] = params.Km_ij
    w.attrib["units"] = "mmol/m3; half-sat"
    
    w = defVar(f, "m_lb", Float64, ("nb",))
    w[:,:] = params.m_lb
    w.attrib["units"] = "m3/mmol; linear death rate of b"

    w = defVar(f, "m_qb", Float64, ("nb",))
    w[:,:] = params.m_qb
    w.attrib["units"] = "m3/mmol; quadratic death rate of b"
    
    w = defVar(f, "K_g", Float64, ("nz",))
    w[:] = params.K_g
    w.attrib["units"] = "m3/mmol; half-sat rate of z"
    
    w = defVar(f, "γ", Float64, ("nz",))
    w[:] = params.γ
    w.attrib["units"] = "fraction of digestion of z"
    
    w = defVar(f, "m_lz", Float64, ("nz",))
    w[:,:] = params.m_lz
    w.attrib["units"] = "m3/mmol; linear death rate of z"

    w = defVar(f, "m_qz", Float64, ("nz",))
    w[:,:] = params.m_qz
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