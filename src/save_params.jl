struct SavePrms
    H::Int64                       # total depth (m)
    dz::Int64                      # height of each box/parcel of water (m)
    np::Int64                      # number of phytoplankton
    nb::Int64                      # number of bacteria
    nz::Int64                      # number of zooplankton
    nn::Int64                      # number of inorganic matter pools
    nd::Int64                      # number of organic matter pools
    pIC::Array{Float64,2}          # initial condition for np p at time 0
    bIC::Array{Float64,2}          # initial condition for nb b at time 0
    zIC::Array{Float64,2}          # initial condition for nz z at time 0
    nIC::Array{Float64,2}          # initial condition for nn n at time 0
    dIC::Array{Float64,2}          # initial condition for nd d at time 0
    oIC::Array{Float64,2}          # initial condition for no o at time 0
    umax_i::Array{Float64,1}       # max uptake rate overall for p_i
    umax_ij::Array{Float64,2}      # max uptake rate of phyto j on n
    Kp_i::Array{Float64,1}         # half saturation rate for p_i 
    Kp_ij::Array{Float64,2}        # half saturation rate of phyto j on n 
    m_lp::Array{Float64,1}         # linear mort of p (1/d) 
    m_qp::Array{Float64,1}         # quadratic mort of p (m3/mmol/d)
    light::Array{Float64, 1}       # light (irradiance) W/m2
    temp_fun::Array{Float64, 1}    # temperature function (modification to metabolic rates)
    K_I::Float64                   # light half-saturation constant
    CMp::Array{Bool,2}             # consumption matrix: nn X np
    vmax_i::Array{Float64,1}       # max uptake rate overall for d_i
    vmax_ij::Array{Float64,2}      # max uptake rate of bacteria j on d_i
    Km_i::Array{Float64,1}         # half saturation rate overall for d_i 
    Km_ij::Array{Float64,2}        # half saturation rate of bacteria j on d_i 
    y_ij::Array{Float64,2}         # yield rate of bacteria j on d_i
    m_lb::Array{Float64,1}         # linear mort of b
    m_qb::Array{Float64,1}         # quadratic mort of b
    prob_generate_d::Array{Float64, 1} # distribution of OM to each d pool from mortality
    CM::Array{Bool,2}              # consumption matrix: nd X nb
    g_max::Array{Float64,1}        # max grazing rate of z
    K_g::Array{Float64,1}          # half saturation rate of z
    γ::Array{Float64,1}            # fraction of assimilation (assimilation efficiency) for z 
    m_lz::Array{Float64,1}         # linear mort of z
    m_qz::Array{Float64,1}         # quadratic mort of z
    GrM::Array{Bool,2}             # grazing matrix: nz x (np + nb)
    pen::Array{Float64,1}          # generalists penalty
    kappa_z::Array{Float64, 1}     # vertical eddy diffusivities
    wd::Array{Float64, 2}          # sinking rate for POM
    ngrid::Int64                   # number of boxes at depth
    e_o::Float64                   # production of O2 (excretion) (mol O2 per mol N uptake)
    yo_ij::Array{Float64, 2}       # O2 yield rate of bacteria j on d_i (mol B/mol O2)
    koverh::Float64                # gas transfer coefficient for each box comprising the mixed layer
    o2_sat::Float64                 # O2 half-sat constant (mmol/m3)
    ml_boxes::Int64                # num boxes in mixed layer, close to 100m total sum
    t_o2relax::Float64             # deep oxygen relaxation (1/day)
    o2_deep::Float64               # mmol/m3
end

function save_params(prms)

    prms_to_save = SavePrms(
                    prms.H, prms.dz, prms.np, prms.nb, prms.nz, prms.nn, prms.nd, prms.pIC, prms.bIC, prms.zIC, prms.nIC, prms.dIC, prms.oIC, 
                    prms.umax_i, prms.umax_ij, prms.Kp_i, prms.Kp_ij, prms.m_lp, prms.m_qp, prms.light, prms.temp_fun, prms.K_I, prms.CMp,
                    prms.vmax_i, prms.vmax_ij, prms.Km_i, prms.Km_ij, prms.y_ij, prms.m_lb, prms.m_qb, prms.prob_generate_d, prms.CM,
                    prms.g_max, prms.K_g, prms.γ, prms.m_lz, prms.m_qz, prms.GrM, prms.pen, prms.kappa_z, prms.wd, prms.ngrid, 
                    prms.e_o, prms.yo_ij, prms.koverh, prms.o2_sat, prms.ml_boxes, prms.t_o2relax, prms.o2_deep
                )

    savefile = replace(prms.fsaven, "out_$(years)y_" => "", ".nc" => "", "results/outfiles/" => "")

    jldopen("results/saved_params/$(savefile).jdl", "w") do f
            write(f, "A", prms_to_save)  
        end

end


function load_saved_params(dt, tt, nrec, nt, fsaven, filename)

    p = jldopen("results/saved_params/$(filename)", "r") do f
                read(f, "A")
    end

    zc = [p.dz/2 : p.dz : p.H - p.dz/2;]  
    zf = [0 : p.dz : p.H;]  

    season = request(message("SE2"), RadioMenu(message("SE1")))
    if season == 1
        mlz = 30
        temp = 6.5 .*exp.(-zc ./ 150) .+ 9 .*exp.(-zc ./ 500) .+ 3
    else
        mlz = 15
        temp = 10 .*exp.(-zc ./ 150) .+ 9 .*exp.(-zc ./ 500) .+ 2.9
    end  
            
    kappa_z = (1e-2  .* exp.(-zf/mlz) .+ 1e-2 .+ 1e-4 .* exp.((zf .- p.H) / 100.)) .* 3600 .* 24 
    kappa_z[1], kappa_z[end] = 0, 0
    coeffs = [0.8, -4000, 293.15, 273.15]
    temp_fun = coeffs[1] .* exp.(coeffs[2] .*(1 ./ (temp .+ coeffs[4]) .- 1 ./ coeffs[3]))

    params = Prms(
        tt, dt, nt, nrec, p.H, p.dz, p.np, p.nb, p.nz, p.nn, p.nd, p.pIC, p.bIC, p.zIC, p.nIC, p.dIC, p.oIC, 
        p.umax_i, p.umax_ij, p.Kp_i, p.Kp_ij, p.m_lp, p.m_qp, p.light, temp_fun, p.K_I, p.CMp,
        p.vmax_i, p.vmax_ij, p.Km_i, p.Km_ij, p.y_ij, p.m_lb, p.m_qb, p.prob_generate_d, p.CM,
        p.g_max, p.K_g, p.γ, p.m_lz, p.m_qz, p.GrM, p.pen, kappa_z, p.wd, p.ngrid, 
        p.e_o, p.yo_ij, p.koverh, p.o2_sat, p.ml_boxes, p.t_o2relax, p.o2_deep, fsaven
    )

    return params, season

end