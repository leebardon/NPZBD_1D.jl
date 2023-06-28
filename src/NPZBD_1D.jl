module NPZBD_1D
    
    import REPL
    using REPL.TerminalMenus
    using Logging, LoggingExtras, Dates
    using SparseArrays, Distributions

    using Distributed
    addprocs(14, exeflags = "--project=$(Base.active_project())")
    println("\n > Number of cores: ", nprocs())
    println(" > Number of workers: ", nworkers())

    logger = TeeLogger(
        MinLevelLogger(FileLogger("logs/info.log"), Logging.Info),
        MinLevelLogger(FileLogger("logs/error.log"), Logging.Warn),
    );

    global_logger(logger)

    include("model.jl")
    include("params.jl")
    include("physics.jl")
    include("utils.jl")
    include("consumption_matrix.jl")
    include("grazing_matrix.jl")
    include("traits.jl")
    include("integrate.jl")
    include("plot.jl")


    #------------------------------------------------------------------------------------------------------------#
    #TODO add test option
    # Figure out how to pause and wait for subroutine to end if test is selected, then quit program, rather than 
    # include("../test/data/prms_for_tests.jl")
    # test_prms = TestPrms()
    # N, P, Z, B, D, track_time, fsaven = run_NPZBD(test_prms, 1)
    #------------------------------------------------------------------------------------------------------------#

    fsave = "results/outfiles/out_1D"
    
    #------------------------------------------------------------------------------------------------------------#
    #   GRID SETUP
    #------------------------------------------------------------------------------------------------------------#
        H = 890                         # depth at SPOT (m)
        dz = 10                         # height per box
        ngrid = Int(H/dz)               # number of boxes
        zc = [dz/2 : dz : H - dz/2;]    # centered depth 
        zf = [0 : dz : H;]              # face depth; 1 longer than zc
    
    #------------------------------------------------------------------------------------------------------------#
    #   COLLECT USER INPUT
    #------------------------------------------------------------------------------------------------------------#
        #TODO expand code to allow for multiple n 
        println(message("START"))
        println(message("DEF"))
        defaults = request(message("DF2"), RadioMenu(message("DF1")))

        if defaults == 1 
            tt, nrec, nd, nb, np, nz, nn, y_i, supply_weight, vmax_i, umax_i, pulse = get_defaults()
        else 
            tt, nrec, nd, nb, np, nz, nn, yield, supply_weight, uptake, uptake_p, pulse = user_select()
            yield == 1 ? y_i = ones(nd)*0.3 : y_i = rand(nd)*0.5
            uptake == 1 ? vmax_i = ordered_uptake_arr(nd) : vmax_i = random_uptake_arr(nd)
            uptake_p == 1 ? umax_i = fill(1., np) : umax_i = fill(1., np) 
            # uptake_p == 1 ? umax_i = ordered_uptake_arr(np, true) : umax_i = random_uptake_arr(np)
        end
            

    # -----------------------------------------------------------------------------------------------------------#
    #   HETEROTROPHIC MICROBE PARAMS
    #------------------------------------------------------------------------------------------------------------#
        if isfile("NPZBD_1D/results/matrices/CM_$(nd)d$(nb)b.jdl")
            CM = load_matrix("CM", nd, nb)
        else
            CM = build_consumption_matrix(nd, nb)
        end

        y_ij = broadcast(*, y_i, CM)
        yo_ij = y_ij*10                     # PLACEHOLDER VALUE mol B/mol O2. not realistic
        num_uptakes = sum(CM, dims=1)[1, :]
        pen = 1 ./ num_uptakes
        Km_i = vmax_i./10 

        println(message("MIC"))

        tradeoff = request(message("TB2"), RadioMenu(message("T1")))  
        if tradeoff == 1 
            vmax_ij, Km_ij = apply_tradeoff(nb, nd, CM, vmax_i)
        else
            vmax_ij = ones(nd, nb) * vmax_i
            Km_ij = ones(nd, nb) * Km_i
        end


    # -----------------------------------------------------------------------------------------------------------#
    #   PHYTOPLANKTON PARAMS 
    #------------------------------------------------------------------------------------------------------------#
        if isfile("NPZBD_1D/results/matrices/CMp_$(np)p$(nd)d$(nb)b$(nz)z.jdl")
            CMp = load_matrix("CMp", nd, nb, nn, np)
        else
            CMp = build_consumption_matrix(nn, np)
        end

        K_I = 10.0         
        e_o = 150/16   
        Kp_i = umax_i./10 

        tradeoff = request(message("TP2"), RadioMenu(message("T1")))  
        if tradeoff == 1 
            umax_ij, Kp_ij = apply_tradeoff(np, nn, CMp, umax_i)
        else
            umax_ij = ones(nn, np) * umax_i
            Kp_ij = ones(nn, np) * Kp_i
        end 

    # -----------------------------------------------------------------------------------------------------------#
    #   ZOOPLANKTON GRAZING 
    #------------------------------------------------------------------------------------------------------------#
        g_max = ones(nz)
        K_g = ones(nz)*0.2
        γ = ones(nz)*0.3

        if isfile("NPZBD_1D/results/matrices/GrM_$(np)p$(nb)b$(nz)z.jdl")
            GrM = load_matrix("GrM", nd, nb, nn, np, nz)
        else
            GrM = build_grazing_matrix(np, nb, nz)
        end

        println(message("CM", nd, nb, nn, np))
        println(message("GM", nd, nb, nn, np, nz))
    # -----------------------------------------------------------------------------------------------------------#
    #   MORTALITY
    #------------------------------------------------------------------------------------------------------------#
        m_lp = ones(np) * 1e-1  
        m_qp = ones(np) * 0.1  # (.1 if grazers, if not, 1)

        m_lb = ones(nb) * 1e-2 
        m_qb = ones(nb) * 0.1 
        # m_qb[1] = 1  # POM consumer 

        m_lz = ones(nz) * 1e-2
        m_qz = ones(nz) * 0.1 

    # -----------------------------------------------------------------------------------------------------------#
    #   ORGANIC MATTER
    #------------------------------------------------------------------------------------------------------------#
        # Distribution of OM from mortality to detritus pools
        if supply_weight == 1 
            prob_generate_d = ones(nd) * (1/nd)
        else
            dist = LogNormal(1.5,2)
            x = rand(dist, nd)
            prob_generate_d =  x / sum(x)
        end

        # Sinking rate for POM  #NOTE could be randomly assigned range 1 t0 10
        ws = zeros(nd)                  # sinking speed of POM (m/day)
        ws[1] = 5
        #NOTE because I have something sinking too fast and accumulating
        # make one pool from average vmax and one could be sinking speed of zero
        # we can remove from ADVECT 
        w = zeros(ngrid + 1)            # water vertical velocity, if there was any
        wd = transpose(repeat(ws, 1, ngrid + 1)) + repeat(w, 1, nd) # ngrid+1 x nd
        wd[1,:] .= 0                    # no flux boundary at surface 
        wd[end,:] .= 0                  # no flux boundary (bottom box accumulates D)

    #------------------------------------------------------------------------------------------------------------#
    #   PHYSICAL ENVIRONMENT
    #------------------------------------------------------------------------------------------------------------#

        # VERTICAL MIXING #NOTE Summer vs winter - mixed layer shallower in summer - try 15 summer, 25 winter (instead of pulse?)
            if pulse == 2 
                mlz = 15     # mixed layer lengthscale
            elseif pulse == 1 || pulse == 3
                mlz = 25
            end

            kappazmin = 1e-4              # min mixing coeff -value for most of the deep ocean (higher at top and bottom)
            kappazmax = 1e-2              # max mixing coeff -value at top of mixed layer (and bottom boundary mixed layer)
            kappa_z = (kappazmax .* exp.(-zf/mlz) .+ kappazmin .+ kappazmax .* exp.((zf .- H) / 100.)) .* 3600 .* 24 
            kappa_z[1] = 0
            kappa_z[end] = 0

        # LIGHT (Irradiance, I) 
            euz = 25                    # euphotic zone lengthscale #NOTE scales with amount of biomass but VERY sensitive
            light_top = 700             # avg incoming PAR = (1400/2)  Light_avg*(cos(t*dt*2*3.1416)+1) for light daily cycle
            light = light_top .* exp.( -zc ./ euz)

        # OXYGEN (air-sea exchange)
            o2_sat = 212.1              # mmol/m3 from calc_oxsat(25+273,35) in matlab. WOCE clim-avg surf T at 10S, E. Pac.
            Kgast = 3e-5*86400          # m/d
            ml_boxes = 100/dz           # discrete n of boxes in the mixed layer, close to 100m total sum
            koverh = Kgast/100/ml_boxes # gas transfer coeff for each of the n boxes comprising the ml. why the 100 here?

        # OXYGEN (deep oxygen relaxation)
            o2_deep = 200.0             # mmol/m3, avg (for ~7 C) and 35
            t_o2relax = 0.01            # 1/day, range from 0.01 to 0.1. Set to 0 to turn off.


        # TEMPERATURE (SPOT along water column)
        #fit to SPOT data (approx 20 to 4, approx 16 to 4)
            if pulse == 2 
                temp = 6.5 .*exp.(-zc ./ 150) .+ 9 .*exp.(-zc ./ 500) .+ 3
            elseif pulse == 1 || pulse == 3
                temp = 10 .*exp.(-zc ./ 150) .+ 9 .*exp.(-zc ./ 500) .+ 2.9
            end

        # TEMPERATURE (modification to metabolic rates)
            temp_coeff_arr = 0.8
            temp_ae_arr = -4000
            temp_ref_arr = 293.15   
            t_kel = 273.15
            temp_fun = temp_coeff_arr .* exp.(temp_ae_arr .*(1 ./ (temp .+ t_kel) .- 1 ./ temp_ref_arr))


    # -----------------------------------------------------------------------------------------------------------#
    #   INITIAL CONDITIONS
    #------------------------------------------------------------------------------------------------------------#
        nIC = ones(Float64, ngrid, nn) * 10.0
        pIC = ones(Float64, ngrid, np) * 0.1 
        zIC = ones(Float64, ngrid, nz) * 0.01
        # dIC = ones(Float64, ngrid, nd) * (0.01 / nd)
        # bIC = ones(Float64, ngrid, nb) * (0.1 / nb)
        dIC = ones(Float64, ngrid, nd) * 0.1
        bIC = ones(Float64, ngrid, nb) * 0.01
        oIC = ones(Float64, ngrid, 1)  * 100.0



    # -----------------------------------------------------------------------------------------------------------#
    #   INSTANTIATE PARAMS & RUN MODEL
    #------------------------------------------------------------------------------------------------------------#

    dt = 0.01
    params = Prms(
                tt, dt, nrec, H, dz, np, nb, nz, nn, nd, pIC, bIC, zIC, nIC, dIC, oIC, 
                umax_i, umax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp,
                vmax_i, vmax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, prob_generate_d, CM,
                g_max, K_g, γ, m_lz, m_qz, GrM, pen, kappa_z, wd, ngrid, 
                e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, fsave
            )
    # params = Prms(
    #             tt, dt, nrec, H, dz, np, nb, nz, nn, nd, pIC, bIC, zIC, nIC, dIC, oIC, 
    #             umax_p, K_n, m_lp, m_qp, CM, y_ij, vmax_i, vmax_ij, Km_i, Km_ij, 
    #             m_lb, m_qb, g_max, K_g, γ, m_lz, m_qz, fsave, GrM, pen, 
    #             prob_generate_d, kappa_z, wd, light, temp_fun, K_I, ngrid, 
    #             e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep,
    #         )

    @info("Model Params: \n $params \n")

    N, P, Z, B, D, O, track_time, fsaven = run_NPZBD(params, pulse)

    save_matrices(CM, CMp, GrM, nd, nb, nn, np, nz)

    depth_plots(fsaven, pulse)

end


#####################################################
#TODO
# """

# 1. Consumption matrix and a function for pulsing nutrients - random and periodic - done

# 2. Grazing and cell size (start with just one grazer for each)

# Add a way to scale the input of P into the consumption matrix, in the way we do for the bacteria
# Later - ammonia oxidizers and nitrite oxidisers - can wait until later in the year 

# """