module NPZBD_1D
    
    import REPL
    using REPL.TerminalMenus
    using Logging, LoggingExtras, Dates
    using SparseArrays, Distributions

    using Distributed
    addprocs(14, exeflags = "--project=$(Base.active_project())")
    println("\n > Number of cores: ", nprocs())
    println(" > Number of workers: ", nworkers())

    include("model.jl")
    include("params.jl")
    include("physics.jl")
    include("utils.jl")
    include("consumption_matrix.jl")
    include("grazing_matrix.jl")
    include("traits.jl")
    include("integrate.jl")
    include("plot.jl")
    include("save_params.jl")


    #------------------------------------------------------------------------------------------------------------#
    #TODO add test option
    # Figure out how to pause and wait for subroutine to end if test is selected, then quit program, rather than 
    # include("../test/data/prms_for_tests.jl")
    # test_prms = TestPrms()
    # N, P, Z, B, D, track_time, fsaven = run_NPZBD(test_prms, 1)
    #------------------------------------------------------------------------------------------------------------#
    


    #------------------------------------------------------------------------------------------------------------#
    #   TIMES AND LOGS
    #------------------------------------------------------------------------------------------------------------#
        println(message("START"))

        simulation_time = request(message("TM2"), RadioMenu(message("TM1")))
        if simulation_time == 1
            years = 1
            tt = 366
            nrec = 7320
        elseif simulation_time == 2
            years = 10
            tt = 3660
            nrec = 73200
        else 
            years = 100
            tt = 36600
            nrec = 732000
        end

        dt = 0.01
        nt = Int(tt/dt)

        fsaven, logger = set_savefiles(now(), years)
    
    #------------------------------------------------------------------------------------------------------------#
    #   COLLECT USER INPUT
    #------------------------------------------------------------------------------------------------------------#
        run_type = request(message("ST2"), RadioMenu(message("ST1")))

        if run_type == 1 
            nd, nb, np, nz, nn, yield, supply_weight, uptake, uptake_p, season = user_select()
            yield == 1 ? y_i = ones(nd)*0.3 : y_i = rand(nd)*0.5
            uptake == 1 ? vmax_i = ordered_uptake_arr(nd) : vmax_i = random_uptake_arr(nd)
            uptake_p == 1 ? umax_i = fill(1., np) : umax_i = random_uptake_arr(np)
        else
            options = readdir("results/saved_params/")
            file = request(message("LVP"), RadioMenu(options))
            params, season = load_saved_params(dt, tt, nrec, nt, fsaven, options[file])

            N, P, Z, B, D, O, track_time = run_NPZBD(params, season)
        
            depth_plots(fsaven, season, years)

            exit()
        end
            
    #------------------------------------------------------------------------------------------------------------#
    #   GRID SETUP
    #------------------------------------------------------------------------------------------------------------#
        H = 890                         # depth at SPOT (m)
        dz = 10                         # height per box
        ngrid = Int(H/dz)               # number of boxes
        zc = [dz/2 : dz : H - dz/2;]    # centered depth 
        zf = [0 : dz : H;]              # face depth; 1 longer than zc

    # -----------------------------------------------------------------------------------------------------------#
    #   HETEROTROPHIC MICROBE PARAMS
    #------------------------------------------------------------------------------------------------------------#
        CM = get_matrix("CM", nd, nb, nn, np, nz)

        y_ij = broadcast(*, y_i, CM)
        yo_ij = y_ij*10                     # PLACEHOLDER VALUE mol B/mol O2. not realistic
        num_uptakes = sum(CM, dims=1)[1, :]
        pen = 1 ./ num_uptakes
        Km_i = vmax_i./10 

        println(message("MIC"))

        tradeoff_b = request(message("TB2"), RadioMenu(message("T1")))  
        if tradeoff_b == 1 
            vmax_ij, Km_ij = apply_tradeoff(nb, nd, CM, vmax_i)
        else
            vmax_ij = ones(nd, nb) * vmax_i
            Km_ij = ones(nd, nb) * Km_i
        end


    # -----------------------------------------------------------------------------------------------------------#
    #   PHYTOPLANKTON PARAMS 
    #------------------------------------------------------------------------------------------------------------#
        CMp = get_matrix("CMp", nd, nb, nn, np, nz)

        K_I = 10.0         
        e_o = 150/16   
        Kp_i = umax_i./10 

        tradeoff_p = request(message("TP2"), RadioMenu(message("T1")))  
        if tradeoff_p == 1 
            umax_ij, Kp_ij = apply_tradeoff(np, nn, CMp, umax_i)
        else
            umax_ij = ones(nn, np) * umax_i
            Kp_ij = ones(nn, np) * Kp_i
        end 

    # -----------------------------------------------------------------------------------------------------------#
    #   ZOOPLANKTON GRAZING 
    #------------------------------------------------------------------------------------------------------------#
        GrM = get_matrix("GrM", nd, nb, nn, np, nz)
    
        g_max = ones(nz)
        K_g = ones(nz)*1.0
        γ = ones(nz)*0.3

    # -----------------------------------------------------------------------------------------------------------#
    #   MORTALITY
    #------------------------------------------------------------------------------------------------------------#
        m_lp = ones(np) * 1e-1  
        m_qp = ones(np) * 0.1  # (.1 if grazers, if not, 1)

        m_lb = ones(nb) * 1e-2 
        m_qb = ones(nb) * 0.1 
        m_qb[1] = 1  # POM consumer 

        m_lz = ones(nz) * 1e-2
        m_qz = ones(nz) * 1.0 

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
        ws[1] = 3 
        w = zeros(ngrid + 1)            # water vertical velocity, if there was any
        wd = transpose(repeat(ws, 1, ngrid + 1)) + repeat(w, 1, nd) # ngrid+1 x nd
        wd[1,:] .= 0                    # no flux boundary at surface 
        wd[end,:] .= 0                  # no flux boundary (bottom box accumulates D)

    #------------------------------------------------------------------------------------------------------------#
    #   PHYSICAL ENVIRONMENT
    #------------------------------------------------------------------------------------------------------------#

        # VERTICAL MIXING 
            season == 1 ? mlz = 30 : mlz = 15  # mixed layer lengthscale
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
            koverh = Kgast/ml_boxes # gas transfer coeff for each of the n boxes comprising the ml. 

        # OXYGEN (deep oxygen relaxation)
            o2_deep = 200.0             # mmol/m3, avg (for ~7 C) and 35
            t_o2relax = 0.01            # 1/day, range from 0.01 to 0.1. Set to 0 to turn off.


        # TEMPERATURE (SPOT along water column)
        #fit to SPOT data (approx 20 to 4, approx 16 to 4)
            if season == 1 
                temp = 6.5 .*exp.(-zc ./ 150) .+ 9 .*exp.(-zc ./ 500) .+ 3
            else
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
        nIC = ones(Float64, ngrid, nn) * 20.0
        pIC = ones(Float64, ngrid, np) * 0.1 
        zIC = ones(Float64, ngrid, nz) * 0.01
        dIC = ones(Float64, ngrid, nd) * (0.01 / nd)
        bIC = ones(Float64, ngrid, nb) * (0.1 / nb)
        oIC = ones(Float64, ngrid, 1)  * 100.0


    # -----------------------------------------------------------------------------------------------------------#
    #   INSTANTIATE PARAMS & RUN MODEL
    #------------------------------------------------------------------------------------------------------------#
        params = Prms(
                    tt, dt, nt, nrec, H, dz, np, nb, nz, nn, nd, pIC, bIC, zIC, nIC, dIC, oIC, 
                    umax_i, umax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp,
                    vmax_i, vmax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, prob_generate_d, CM,
                    g_max, K_g, γ, m_lz, m_qz, GrM, pen, kappa_z, wd, ngrid, 
                    e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, fsaven
                )

        @info("Model Params: \n $params \n") ; print_info(params)

        N, P, Z, B, D, O, track_time = run_NPZBD(params, season)
    
        save_matrices(CM, CMp, GrM, nd, nb, nn, np, nz)
        
        depth_plots(fsaven, season, years)

        save_prm = request(message("SVP2"), RadioMenu(message("SVP1")))
        save_prm == 1 ? save_params(params) : exit()

end


#####################################################
#TODO
# """

# 1. Consumption matrix and a function for pulsing nutrients - random and periodic - done

# 2. Grazing and cell size (start with just one grazer for each)

# Add a way to scale the input of P into the consumption matrix, in the way we do for the bacteria
# Later - ammonia oxidizers and nitrite oxidisers - can wait until later in the year 

# """