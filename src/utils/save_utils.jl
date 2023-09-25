
using DataFrames, NCDatasets, JLD


function set_savefiles(launch_time, season, years, np, nz, nb, nd)

    season == 1 ? season_str = "Wi" : season_str = "Su"
    fsave = "results/outfiles/$(season_str)$(years)y"
    fsaven = string(fsave, "_", Dates.format(launch_time, "yymmdd_HH:MM"), "_$(np)P$(nz)Z$(nb)B$(nd)D.nc")

    return fsaven

end


function define_dims(ds, prms, nrec1)

    vars = Dict("nb" => prms.nb, "nd" => prms.nd, "nrec" => nrec1)

    for (k, v) in x
        defDim(ds, k, v)
    end

    return ds

end


function save_full_run(p, b, z, n, d, o, timet, v, uptake, tst, tfn, prms, season_num)

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
    f.attrib["title"] = "NPZBD 1D model i/o"
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

    # w = defVar(f, "temp_fun", Float64, ("ndepth",))
    # w[:] = prms.temp_fun
    # w.attrib["units"] = "temp mod to metabolic rate"

    # --------------------------------------------------

    w = defVar(f,"uptake",Float64,("nd","nb"))
    w[:,:,:] = uptake
    w.attrib["units"] = "mmol/m3 C per d; uptake matrix"

    w = defVar(f,"SW",Float64,("nd",))
    w[:,:] = prms.prob_generate_d 
    w.attrib["units"] = "Ind C supply weight: probability"
    
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



function save_endpoints(n, p, z, b, d, o, prms, season)

    outdir = "/home/lee/Dropbox/Development/NPZBD_1D/"
    ep_path = replace(prms.fsaven, "results/outfiles" => "results/outfiles/endpoints", ".nc" => "_ep.nc")
    path = joinpath(outdir, ep_path) 
    println("\nSaving endpoints to: ", path)

    ep = get_endpoints([n, p, z, b, d, o])
    n, p, z, b, d, o = ep[1], ep[2], ep[3], ep[4], ep[5], ep[6]

    f = NCDataset(path, "c") 

    # define the dim of p, b, z, n, d
    defDim(f,"np", prms.np)
    defDim(f,"nb",prms.nb)
    defDim(f,"nz",prms.nz)
    defDim(f,"nn",prms.nn)
    defDim(f,"nd",prms.nd)
    defDim(f,"no",1)

    # define the dim of the depth
    defDim(f,"ndepth",89)
    defDim(f,"ndepth1",90)

    # define the dim of the time length
    nrec1 = Int(prms.nrec+1) #bc i added time 0
    nprey = prms.np + prms.nb
    
    defDim(f,"nrec",nrec1)
    defDim(f,"nprey",nprey)
   
    # info
    f.attrib["title"] = "NPZBD 1D model endpoints"
    f.attrib["Season"] = season

    # simulated results
    w = defVar(f,"p",Float64,("ndepth" ,"np"))
    w[:,:] = p
    w.attrib["units"] = "mmol/m3 N biomass"

    w = defVar(f,"b",Float64,("ndepth" ,"nb"))
    w[:,:] = b
    w.attrib["units"] = "mmol/m3 N biomass"

    w = defVar(f,"z",Float64,("ndepth" ,"nz"))
    w[:,:] = z
    w.attrib["units"] = "mmol/m3 N biomass"
    
    w = defVar(f,"n",Float64,("ndepth" ,"nn"))
    w[:,:] = n
    w.attrib["units"] = "mmol/m3 N OM"

    w = defVar(f,"d",Float64,("ndepth" ,"nd"))
    w[:,:] = d
    w.attrib["units"] = "mmol/m3 N OM"
    
    w = defVar(f,"o",Float64,("ndepth" ,"no"))
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

    w = defVar(f, "ws_POM", Float64, ())
    w[:] = prms.ws_POM
    w.attrib["units"] = "sinking rate of POM in m/day"

    w = defVar(f, "temp_fun", Float64, ("ndepth",))
    w[:] = prms.temp_fun
    w.attrib["units"] = "temp mod to metabolic rate"

    # --------------------------------------------------

    # w = defVar(f,"uptake",Float64,("nd","nb"))
    # w[:,:,:] = uptake
    # w.attrib["units"] = "mmol/m3 C per d; uptake matrix"
    
    w = defVar(f, "umax_i", Float64, ("np",))
    w[:] = prms.umax_i
    w.attrib["units"] = "m3/mmol/d; max growth rate of p"

    w = defVar(f, "Fg_p", Float64, ("np",))
    w[:,:] = prms.Fg_p
    w.attrib["units"] = "per p; fraction proteome assigned to growth"

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

    w = defVar(f, "Fg_b", Float64, ("nb",))
    w[:,:] = prms.Fg_b
    w.attrib["units"] = "per b; fraction proteome devoted to growth"
    
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


# ds = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_8P6Z13B5D_230915_00:03.nc")
# save_endpoints(ds,"summer", 100, 8, 13, 6, 1, 5, 1, "Wi100y_8P6Z13B5D_230915_00:03_ep.nc")