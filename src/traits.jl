###################################################
#SUBSTRATE TRAITS ("averages")

function ordered_uptake_arr(n, phy=false)

    max_uptake_i = collect(10 .^ range(-2,stop=1,length=n))

    if phy==true
        max_uptake_i = collect(10 .^ range(0,stop=1,length=n))
    end

    return max_uptake_i

end


function random_uptake_arr(n)
    # randomly selecting along a log range 
    xhi = log10(1e1)
    xlo = log10(1e-1) 
    max_uptake_i = 10 .^ ((rand(n).-0.5).*(xhi-xlo).+mean([xhi xlo]))

    return max_uptake_i

end


function apply_tradeoff(consumers, resources, CM, max_i)
    #fraction of the proteome devoted to growth (F_g) vs. affinity (F_a)
    #NOTE could also explore growth-defense tradeoff https://www.nature.com/articles/s41396-020-0619-1 
    F_g = rand(consumers)
    F_a = 1. .- F_g

    if resources == 1
        max_ij = set_umax_ij(resources, consumers, max_i, CM, F_g)
        K_ij = set_Kp_ij(resources, consumers, F_a, CM, max_ij)
    else 
        max_ij = set_vmax_ij(resources, consumers, max_i, CM, F_g)
        K_ij = set_Km_ij(resources, consumers, F_a, CM, max_ij)
    end 

    return max_ij, K_ij

end


function set_vmax_ij(nd, nb, vmax_i, CM, F_g)
    # if F_g == 0.5 and F_a is 0.5, then max_ij == max_i
    vmax_ij = zeros(nd, nb)
    for i = 1:nd
        vmax_ij[i,:] = (vmax_i[i] * F_g) ./ (0.5 .* CM[i,:])  
    end

    return vmax_ij
end 


function set_umax_ij(nn, np, umax_i, CM, F_g)
    # if F_g == 0.5 and F_a is 0.5, then max_ij == max_i
    umax_ij = zeros(nn, np)
    for i = 1:nn
        umax_ij[i,:] = (umax_i[i] * F_g) ./ (0.5 .* CM[i,:])  
    end

    return umax_ij
end 


function set_Km_ij(nd, nb, F_a, CM, vmax_ij)
    # if F_a is 0.5, then Km_ij == Km_i
    affinity = zeros(nd, nb)
    Km_ij = zeros(nd, nb)
        for i = 1:nd
        affinity[i,:] = 10*F_a ./ 0.5.*CM[i,:]
        Km_ij[i,:] = vmax_ij[i,:] ./ affinity[i,:] 
    end

    return Km_ij
end


function set_Kp_ij(nn, np, F_a, CM, umax_ij)
    # if F_a is 0.5, then Km_ij == Km_i
    affinity = zeros(nn, np)
    Kp_ij = zeros(nn, np)
    for i = 1:nn
        affinity[i,:] = 10*F_a ./ 0.5.*CM[i,:]
        Kp_ij[i,:] = umax_ij[i,:] ./ affinity[i,:] 
    end

    return Kp_ij
end



