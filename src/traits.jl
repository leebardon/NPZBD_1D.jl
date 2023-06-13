###################################################
#SUBSTRATE TRAITS ("averages")

function ordered_vmax(nd)

    vmax_i = collect(10 .^ range(-2,stop=1,length=nd))

    return vmax_i

end


function random_vmax(nd)
    # randomly selecting along a log range 
    xhi = log10(1e1)
    xlo = log10(1e-1) 
    vmax_i = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))

    return vmax_i

end


function growth_affinity_tradeoff(nb, nd, CM, vmax_i)
    #fraction of the proteome devoted to growth (F_g) vs. affinity (F_a)
    F_g = rand(nb)
    F_a = 1. .- F_g

    vmax_ij = set_vmax_ij(nd, nb, vmax_i, CM, F_g)
    Km_ij = set_km_affinity(nd, nb, F_a, CM, vmax_ij)

    return vmax_ij, Km_ij

end


function set_vmax_ij(nd, nb, vmax_i, CM, F_g)
    # if F_g == 0.5 and F_a is 0.5, then vmax_ij == vmax_i
    vmax_ij = zeros(nd,nb)
    for i = 1:nd
        vmax_ij[i,:] = vmax_i[i]*F_g./0.5.*CM[i,:]  
    end

    return vmax_ij

end 


function set_km_affinity(nd, nb, F_a, CM, vmax_ij)
    # if F_a is 0.5, then Km_ij == Km_i
    aff = zeros(nd,nb)
    Km_ij = zeros(nd,nb)
    for i = 1:nd
        aff[i,:] = 10*F_a./0.5.*CM[i,:]
        Km_ij[i,:] = vmax_ij[i,:]./aff[i,:] 
    end

    return Km_ij
end



