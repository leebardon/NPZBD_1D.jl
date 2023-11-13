
include("utils/utils.jl")
include("utils/save_utils.jl")

function run_prescribed_darwin(years, days, nrec, dt, nt, run_type)

    nd = 8
    nb = 13
    np = 6
    nz = 3
    nn = 1

    println(message("ENV"))
    pulse = request(message("P2"), RadioMenu(message("P1")))
    season = request(message("SE2"), RadioMenu(message("SE1")))

    global fsaven = set_savefiles(now(), season, years, np, nz, nb, nd)

    vmax_i = [1.8, 3.0, 5.0, 7.0, 9.0, 12.0]
    vmax_ij = reshape([0.9, 2.1, 5.0, 8.4, 12.24, 18.96], 1, 6)
    Fg_p = [0.1, 0.25, 0.5, 0.68, 0.79, 0.91]

    umax_i = [1., 1., 1., 1., 1., 1., 1., 1.]
    umax_ij =  [5.17  0  0  0  0  0  0  0  0  0  0  0  0    # first 3 are POM, next 10 are DOM
                0  0.92  0  0  0  0  0  0  0  0  0  0  0    
                0  0  0.16  0  0  0  0  0  0  0  0  0  0  
                0  0  0  29.1  0  0  0  0  16.3  0  0  0  0 
                0  0  0  0  5.17  0  0  0  0  2.9  0  0  0 
                0  0  0  0  0  0.92  0  0  0  0  0.51  0  0 
                0  0  0  0  0  0  0.16  0  0  0  0  0.091  0 
                0  0  0  0  0  0  0  0.029  0  0  0  0  0.016 ] 

    Km_ij = [0.28  0  0  0  0  0  0  0  0  0  0  0  0    # in units of K_DON, calculated from K_DOC in zakem paper by dividing by 5
             0  0.05  0  0  0  0  0  0  0  0  0  0  0    
             0  0  0.0088  0  0  0  0  0  0  0  0  0  0  
             0  0  0  1.56  0  0  0  0  0.72  0  0  0  0 
             0  0  0  0  0.28  0  0  0  0  0.128  0  0  0 
             0  0  0  0  0  0.05  0  0  0  0  0.022  0  0 
             0  0  0  0  0  0  0.0088  0  0  0  0  0.004  0 
             0  0  0  0  0  0  0  0.00156  0  0  0  0  0.00072 ]  
    Fg_b = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  #dummy as Fg vals used to set umax_ij/Km_ij, and those vals already provided

    CM = [1  0  0  0  0  0  0  0  0  0  0  0  0   
          0  1  0  0  0  0  0  0  0  0  0  0  0    
          0  0  1  0  0  0  0  0  0  0  0  0  0  
          0  0  0  1  0  0  0  0  1  0  0  0  0 
          0  0  0  0  1  0  0  0  0  1  0  0  0 
          0  0  0  0  0  1  0  0  0  0  1  0  0 
          0  0  0  0  0  0  1  0  0  0  0  1  0 
          0  0  0  0  0  0  0  1  0  0  0  0  1 ] 

    # GrM - first 6 cols are phyto, next 10 are dom consuming bacteria, rows are zoo
    # Z1 - pom consumers; Z2 - larger (copio, 0.6um) and phyto grazer; Z3 - smaller (oligo, 0.3um) 
    GrM = [1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 
           0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0 
           0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0 ] 

    CMp = get_matrix("CMp", nd, nb, nn, np, nz)


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
    y_i = ones(nd)*0.3
    y_ij = broadcast(*, y_i, CM) 
    yo_ij = y_ij*10                     # PLACEHOLDER VALUE mol B/mol O2. not realistic
    num_uptakes = sum(CM, dims=1)[1, :]
    pen = 1 ./ num_uptakes
    Km_i = umax_i./10 


    # -----------------------------------------------------------------------------------------------------------#
    #   PHYTOPLANKTON PARAMS 
    #------------------------------------------------------------------------------------------------------------
    K_I = 10.0         
    e_o = 150/16   
    Kp_i = vmax_i./10 

    Fa_p = 1. .- Fg_p
    vmax_ij = set_vmax_ij(nn, np, vmax_i, Fg_p)
    Kp_ij = set_Kp_ij(nn, np, Fa_p, CMp, vmax_ij)


    # -----------------------------------------------------------------------------------------------------------#
    #   ZOOPLANKTON GRAZING 
    #------------------------------------------------------------------------------------------------------------
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

    m_lz = ones(nz) * 1e-2
    m_qz = ones(nz) * 1.0 


    # -----------------------------------------------------------------------------------------------------------#
    #   ORGANIC MATTER
    #------------------------------------------------------------------------------------------------------------#
    # prob_generate_d calculated as: dist = Lognormal(1.5, 2) 
    # pom = rand(dist, 3) ; dom = rand(dist, 5) ; d = vcat(pom, dom) ; out = d / sum(d)
    prob_generate_d = [0.05, 0.2, 0.05, 0.05, 0.1, 0.40, 0.1, 0.05]

    # Sinking rate for POM  
    ws = zeros(nd)                  
    ws[1], ws[2], ws[3] = 10.0, 5.0, 1.0


    #------------------------------------------------------------------------------------------------------------#
    #   PHYSICAL ENVIRONMENT
    #------------------------------------------------------------------------------------------------------------#

    # VERTICAL VELOCITY
        w = zeros(ngrid + 1)           
        wd = transpose(repeat(ws, 1, ngrid + 1)) + repeat(w, 1, nd) # ngrid+1 x nd
        wd[1,:] .= 0                    # no flux boundary at surface 
        wd[end,:] .= 0                  # no flux boundary (bottom box accumulates D)

    # VERTICAL MIXING 
        season == 1 ? mlz = 25 : mlz = 15 # mixed layer lengthscale
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
            temp = 6.5 .* exp.(-zc ./ 150) .+ 9 .* exp.(-zc ./ 500) .+ 3
        else
            temp = 10 .* exp.(-zc ./ 150) .+ 9 .* exp.(-zc ./ 500) .+ 2.9
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
                g_max, K_g, γ, m_lz, m_qz, GrM, pen, kappa_z, wd, ngrid, pulse, 
                e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, fsaven
            )

    @info("Model Params: \n $params \n") ; print_info(params)

    N, P, Z, B, D, O, track_time = run_NPZBD(params, season)

    # save_matrices(CM, CMp, GrM, nd, nb, nn, np, nz)
    plot_state_vars(fsaven, season)
    # plot_time_series(fsaven, season)
    # rstar_analysis(fsaven, season)
    # copio_index_analysis()
    # diversity_analysis()

    exit()
    
end

