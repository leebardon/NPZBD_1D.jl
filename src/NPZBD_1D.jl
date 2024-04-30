 module NPZBD_1D
    
    import REPL
    using REPL.TerminalMenus
    using Logging, LoggingExtras 
    using Dates, Printf, Parameters
    using SparseArrays, Distributions, LinearAlgebra
    using Statistics, StatsBase, Random, StableRNGs
    using DataFrames, NCDatasets, CSV
    using Plots, ColorSchemes, Colors, LaTeXStrings
    using Distributed

    GlobalRNG = StableRNG(123)
    addprocs(15, exeflags = "--project=$(Base.active_project())")
    println("\n > Number of cores: ", nprocs())
    println(" > Number of workers: ", nworkers())

    include("model.jl")
    include("params.jl")
    include("physics.jl")
    include("consumption_matrix.jl")
    include("grazing_matrix.jl")
    include("traits.jl")
    include("integrate.jl")
    include("prescribed.jl")
    include("prescribed_darwin.jl")
    include("rstar.jl")
    include("copio_indexing.jl")
    include("nutrient_pulse.jl")
    # include("plotting/heatmaps.jl")
    include("plotting/state_var_plots.jl")
    include("plotting/timeseries_plots.jl")
    include("plotting/rstar_plots.jl")
    include("plotting/state_var_dar.jl")
    include("utils/save_params.jl")
    include("utils/utils.jl")
    include("utils/save_utils.jl")

    # include("plot_test.jl")

    #TODO write tests to check all are equal?
    #------------------------------------------------------------------------------------------------------------#
    #   TIMES AND LOGS
    #------------------------------------------------------------------------------------------------------------#
        println(message("START"))

        run_type = request(message("ST2"), RadioMenu(message("ST1")))

        simulation_time = request(message("TM2"), RadioMenu(message("TM1")))
        if simulation_time == 1
            years = 0
            days = 365
            nrec = 73000
        elseif simulation_time == 2
            years = 2
            days = 732
            nrec = 14640
        elseif simulation_time == 3
            years = 10
            days = 3660
            nrec = 73200
        elseif simulation_time == 4
            years = 30
            days = 10980
            nrec = 219600
        elseif simulation_time == 5
            years = 50
            days = 18300
            nrec = 366000
        else 
            years = 100
            days = 36600
            nrec = 732000
        end

        # nrec = 20 per day when dt = 0.01 (i.e. 100 ts per day, tracking recorded every 5 ts)
        simulation_time == 1 ? dt = 0.001 : dt = 0.01
        nt = Int(days/dt)

        logger = set_logger(now())
    
    #------------------------------------------------------------------------------------------------------------#
    #   COLLECT USER INPUT
    #------------------------------------------------------------------------------------------------------------#

        if run_type == 1 
            # Randomly generated params and consumption/grazing matrices according to user selected state var ratios
            nd, nb, np, nz, nn, y_i, supply_weight, umax_i, vmax_i, season, pulse = user_select(run_type)
            CM = get_matrix("CM", nd, nb, nn, np, nz)
            GrM = get_matrix("GrM", nd, nb, nn, np, nz)
            CMp = get_matrix("CMp", nd, nb, nn, np, nz)

        elseif run_type == 2
            darwin = request(message("DR2"), RadioMenu(message("DR1")))
            darwin == 1 ? run_prescribed_darwin(years, days, nrec, dt, nt, run_type) : nothing
            
            # Loads prescribed params and matrices from prescribed.jl
            nd, nb, np, nz, nn, y_i, supply_weight, umax_i, vmax_i, season, pulse = user_select(run_type)
            CM = get_prescribed_params("CM") 
            GrM = get_prescribed_params("GrM") 
            CMp = get_matrix("CMp", nd, nb, nn, np, nz)

        elseif run_type == 3
            # Continues from end of previous runs (prompts user to select .nc file in results/outfiles )
            continue_type = request(message("CT2"), RadioMenu(message("CT1")))

            H, dz, nIC, pIC, zIC, bIC, dIC, oIC, nn, np, nz, nb, nd, 
            vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p, 
            umax_i, umax_ij, Km_i, Km_ij, m_lb, m_qb, y_ij, prob_generate_d, CM, Fg_b,
            g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, e_o, yo_ij,
            koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, season, prev_fname = get_previous_params(continue_type)

            fsaven = continuation_savefile(prev_fname)
            params = Prms(
                        days, dt, nt, nrec, H, dz, np, nb, nz, nn, nd, pIC, bIC, zIC, nIC, dIC, oIC, 
                        vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p,
                        umax_i, umax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, prob_generate_d, CM, Fg_b,
                        g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
                        e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, fsaven
            )

            log_params(params, season)
            N, P, Z, B, D, O, track_time = run_NPZBD(params, season)
            exit()

        elseif run_type == 4
            # Takes output of previous run and tracks single pulse over 365 days
            H, dz, nIC, pIC, zIC, bIC, dIC, oIC, nn, np, nz, nb, nd, 
            vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p, 
            umax_i, umax_ij, Km_i, Km_ij, m_lb, m_qb, y_ij, prob_generate_d, CM, Fg_b,
            g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, e_o, yo_ij,
            koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, season, prev_fname = get_previous_params()

            bIC[:,4:6,:] = bIC[:,9:11,:]
            bIC[:,12,:] = bIC[:,7,:]
            bIC[:,13,:] = bIC[:,8,:]

            # dIC[1:5,:,:] .= 1.0

            years = 0
            bloom=true
            fsaven = continuation_savefile(prev_fname, bloom)
            params = Prms(
                        days, dt, nt, nrec, H, dz, np, nb, nz, nn, nd, pIC, bIC, zIC, nIC, dIC, oIC, 
                        vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p,
                        umax_i, umax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, prob_generate_d, CM, Fg_b,
                        g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
                        e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, fsaven
            )

            log_params(params, season)
            N, P, Z, B, D, O, track_time = run_NPZBD(params, season, bloom)     
            exit()

        end

        fsaven = set_savefiles(now(), season, pulse, years, np, nz, nb, nd)
        
        
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
        run_type != 3 ? y_ij = broadcast(*, y_i, CM) : nothing
        yo_ij = y_ij*10                     # PLACEHOLDER VALUE mol B/mol O2. not realistic
        num_uptakes = sum(CM, dims=1)[1, :]
        pen = 1 ./ num_uptakes
        Km_i = umax_i./10 

        if run_type != 3
            println(message("MIC"))
            tradeoff_b = request(message("TB2"), RadioMenu(message("T1")))  
            if tradeoff_b == 1 
                umax_ij, Km_ij, Fg_b = apply_tradeoff(nb, nd, CM, umax_i, season, run_type)
            else
                umax_ij = ones(nd, nb) * umax_i
                Km_ij = ones(nd, nb) * Km_i
            end
        else; end


    # -----------------------------------------------------------------------------------------------------------#
    #   PHYTOPLANKTON PARAMS 
    #------------------------------------------------------------------------------------------------------------
        Kp_i = vmax_i./10 

        if run_type != 3
            tradeoff_p = request(message("TP2"), RadioMenu(message("T1")))  
            if tradeoff_p == 1 
                vmax_ij, Kp_ij, Fg_p = apply_tradeoff(np, nn, CMp, vmax_i, season, run_type)
            else
                vmax_ij = ones(nn, np) * vmax_i
                Kp_ij = ones(nn, np) * Kp_i
            end 
        else; end

        # save_prm = request(message("SVP2"), RadioMenu(message("SVP1")))


    # -----------------------------------------------------------------------------------------------------------#
    #   ZOOPLANKTON GRAZING 
    #------------------------------------------------------------------------------------------------------------
        g_max = ones(nz)
        K_g = ones(nz)*1.0
        γ = ones(nz)*0.3


    # -----------------------------------------------------------------------------------------------------------#
    #   MORTALITY RATES (m3/mmol/day)
    #------------------------------------------------------------------------------------------------------------#
        m_lp = ones(np) * 0.1 
        m_qp = ones(np) * 0.1  # (.1 if grazers, if not, 1)

        m_lb = ones(nb) * 0.01 
        m_qb = ones(nb) * 0.1 
        # m_qb[1] = 1  # POM consumer (use if no pom grazer)

        m_lz = ones(nz) * 0.01
        m_qz = ones(nz) * 1.0 


    # -----------------------------------------------------------------------------------------------------------#
    #   ORGANIC MATTER
    #------------------------------------------------------------------------------------------------------------#
        # Distribution of OM from mortality to detritus pools
        if run_type != 3
            if supply_weight == 1 
                prob_generate_d = get_prescribed_params("supply_weight") 
            elseif supply_weight == 2
                prob_generate_d = ones(nd) * (1/nd)
            else 
                dist = LogNormal(1.5,2)
                x = rand(dist, nd)
                prob_generate_d =  x / sum(x)
            end
        else; end

        # Sinking rate for POM  
        ws = zeros(nd)                  
        #NOTE Uncomment for 3 POM (m/day)
        ws_POM1, ws_POM2, ws_POM3 = 6.0, 8.0, 10.0
        ws[1], ws[2], ws[3] = ws_POM1, ws_POM2, ws_POM3

        w = zeros(ngrid + 1)            # water vertical velocity
        wd = transpose(repeat(ws, 1, ngrid + 1)) + repeat(w, 1, nd) # ngrid+1 x nd
        wd[1,:] .= 0                    # no flux boundary at surface 
        wd[end,:] .= 0                  # no flux boundary (bottom box accumulates D)


    #------------------------------------------------------------------------------------------------------------#
    #   PHYSICAL ENVIRONMENT
    #------------------------------------------------------------------------------------------------------------#

        # VERTICAL MIXING 
            season == 1 ? mlz = 25 : mlz = 15 # mixed layer lengthscale
            kappazmin = 1e-4              # min mixing coeff -value for most of the deep ocean (higher at top and bottom)
            kappazmax = 1e-2              # max mixing coeff -value at top of mixed layer (and bottom boundary mixed layer)
            kappa_z = (kappazmax .* exp.(-zf/mlz) .+ kappazmin .+ kappazmax .* exp.((zf .- H) / 100.)) .* 3600 .* 24 
            kappa_z[1] = 0
            kappa_z[end] = 0

        # LIGHT (Irradiance, I) 
            K_I = 10                    # Light half-saturation constant
            euz = 25                    # euphotic zone lengthscale 
            light_top = 700             # avg incoming PAR = (1400/2)  Light_avg*(cos(t*dt*2*3.1416)+1) for light daily cycle
            light = light_top .* exp.( -zc ./ euz)

        # OXYGEN (air-sea exchange)
            o2_sat = 212.1              # mmol/m3 from calc_oxsat(25+273,35) in matlab. WOCE clim-avg surf T at 10S, E. Pac.
            Kgast = 3e-5*86400          # m/d
            ml_boxes = 100/dz           # discrete n of boxes in the mixed layer, close to 100m total sum
            koverh = Kgast/ml_boxes # gas transfer coeff for each of the n boxes comprising the ml. 

        # OXYGEN (deep oxygen relaxation)
            o2_deep = 200.0             # mmol/m3, avg (for ~7 C) and 35
            t_o2relax = 0.01            # 1/day, range from 0.01 to 0.1. Set to 0 to turn off

        # OXYGEN (biotic production)
            yo_ij = y_ij*10             #NOTE PLACEHOLDER VALUE (mol B/mol O2) not realistic
            e_o = 150/16                # production of O2 (excretion). mol O2/mol N uptake


        #NOTE do proper temp fit. Win colder at surface and 5m, spring colder at DCM and below
        # Similarly, summer warmer at surface and 5m, fall warmer at dcm and below
        # TEMPERATURE (SPOT along water column)
            # if season == 1 
            #     temp = 4.2 .* exp.(-zc ./ 150) .+ 8.2 .* exp.(-zc ./ 500) .+ 2.9
            # elseif season == 2
            #     temp = 5.2 .* exp.(-zc ./ 150) .+ 7.5 .* exp.(-zc ./ 500) .+ 3
            # elseif season == 3
            #     temp = 9 .* exp.(-zc ./ 150) .+ 8 .* exp.(-zc ./ 500) .+ 2.9
            # else
            #     temp = 7 .* exp.(-zc ./ 150) .+ 7.5 .* exp.(-zc ./ 500) .+ 3
            # end
            if season == 1 
                temp = 4.2 .* exp.(-zc ./ 150) .+ 8.2 .* exp.(-zc ./ 500) .+ 2.9
            else
                temp = 9 .* exp.(-zc ./ 150) .+ 8 .* exp.(-zc ./ 500) .+ 2.9
            end

        # TEMPERATURE (modification to metabolic rates)
            temp_coeff_arr = 0.8
            temp_ae_arr = -4000
            temp_ref_arr = 293.15   
            t_kel = 273.15
            temp_fun = temp_coeff_arr .* exp.(temp_ae_arr .*(1 ./ (temp .+ t_kel) .- 1 ./ temp_ref_arr))


    # -----------------------------------------------------------------------------------------------------------#
    #   INITIAL CONDITIONS (mmol N/m3)
    #------------------------------------------------------------------------------------------------------------#
        
        if run_type != 3
            nIC = ones(Float64, ngrid, nn) * 20.0
            pIC = ones(Float64, ngrid, np) * 0.1 
            zIC = ones(Float64, ngrid, nz) * 0.01
            dIC = ones(Float64, ngrid, nd) * 0.01
            bIC = ones(Float64, ngrid, nb) * 0.1 
            oIC = ones(Float64, ngrid, 1)  * 100.0
        else
            nIC = n_cont
            pIC = p_cont
            zIC = z_cont
            bIC = b_cont
            dIC = d_cont
            oIC = o_cont
        end


    # -----------------------------------------------------------------------------------------------------------#
    #   INSTANTIATE PARAMS & RUN MODEL
    #------------------------------------------------------------------------------------------------------------#
        params = Prms(
                    days, dt, nt, nrec, H, dz, np, nb, nz, nn, nd, pIC, bIC, zIC, nIC, dIC, oIC, 
                    vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p,
                    umax_i, umax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, prob_generate_d, CM, Fg_b,
                    g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
                    e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, fsaven
                )

        # @info("Model Params: \n $params \n") ; print_info(params)
        log_params(params, season)

        N, P, Z, B, D, O, track_time = run_NPZBD(params, season)
    
        # save_matrices(CM, CMp, GrM, nd, nb, nn, np, nz)
        # plot_state_vars(fsaven, season)
        # plot_time_series(fsaven, season)
        # rstar_analysis(fsaven, season)
        # copio_index_analysis()
        # diversity_analysis()

        # save_prm == 1 ? save_params(params) : exit()
        exit()

end

