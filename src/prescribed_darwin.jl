

function get_prescribed_darwin(param)
    #TODO adapt main script to read this file and eith create new params generator or work out missing stuff 
    # (e.g. work back from umax_ij to get umax_i, same for Km_ij)
    ## ----------------------------------------------------------------------------------------------------------- 
    ##                DARWIN PARAMS ----- 1N, 6P, 3Z, 13B, 8D 
    ## ----------------------------------------------------------------------------------------------------------- 
    # vmax_i = [1.8, 3.0, 5.0, 7.0, 9.0, 12.0]
    if param == "vmax_ij"
        out = [0.9, 2.1, 5.0, 8.4, 12.24, 18.96]

    elseif param == "Kp_ij"
        out = [0.06, 0.161538, 0.5, 1.05, 1.9125, 4.51429]

    elseif param == "umax_ij"
        # first 3 are POM, next 10 are DOM
        out = [5.17, 0.92, 0.16, 29.1, 5.17, 0.92, 0.16, 0.029, 16.3, 2.9, 0.51, 0.091, 0.016]

    elseif param == "Km_ij"
        # in units of K_DON, calculated from K_DOC in zakem paper by dividing by 5
        out = [0.28, 0.05, 0.0088, 1.56, 0.28, 0.05, 0.0088, 0.00156, 0.72, 0.128, 0.022, 0.004, 0.00072]

    elseif param == "supply_weight"
        # calculated as: dist = Lognormal(1.5, 2) 
        # pom = rand(dist, 3) ; dom = rand(dist, 5) ; d = vcat(pom, dom) ; out = d / sum(d)
        # arranged so that pom and dom follow similar pattern as zakem paper
        out = [0.0414275, 0.110539, 0.0522625, 0.00254494, 0.045777, 0.704007, 0.026117, 0.0173254]

    elseif param == "CM"
        # cols are bacteria, rows are OM pools, first 3 rows are POM
        out = [1  0  0  0  0  0  0  0  0  0  0  0  0   
               0  1  0  0  0  0  0  0  0  0  0  0  0    
               0  0  1  0  0  0  0  0  0  0  0  0  0  
               0  0  0  1  0  0  0  0  1  0  0  0  0 
               0  0  0  0  1  0  0  0  0  1  0  0  0 
               0  0  0  0  0  1  0  0  0  0  1  0  0 
               0  0  0  0  0  0  1  0  0  0  0  1  0 
               0  0  0  0  0  0  0  1  0  0  0  0  1 ] 

    elseif param == "GrM"
        # first 6 cols are phyto, next 10 are dom consuming bacteria, rows are zoo
        # Z1 - pom consumers (for simplicity in my current model - zakem uses implicit grazer for pom)
        # Z2 - larger (copio, 0.6um) and phyto grazer
        # Z3 - smaller (oligo, 0.3um) 
        out = [1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 
               0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0 
               0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0 ] 
    
    end
    
    return out

end
