
function advection(c, wd, dz)
    c = vcat(zeros(2,size(c)[2]), c, zeros(2, size(c)[2]))
    positive_wd = (wd + abs.(wd))./2. 
    negative_wd = (wd - abs.(wd))./2.
    f_up = positive_wd[2:end,:].*c[3:end-2,:] + negative_wd[2:end,:].*c[4:end-1,:]
    f_down = positive_wd[1:end-1,:].*c[2:end-3,:] + negative_wd[1:end-1,:].*c[3:end-2,:]
    adv = (f_up .- f_down)./dz 
    return adv 
end


function diffusion(c, kappa_z, dz)
    c = vcat(transpose(zeros(size(c)[2])), c, transpose(zeros(size(c)[2])))
    kappa_z_len = size(kappa_z)[1]
    f_up = kappa_z[1:kappa_z_len-1] .* (c[2:kappa_z_len,:] .- c[1:kappa_z_len-1,:]) ./dz 
    f_down = kappa_z[2:kappa_z_len] .* (c[3:kappa_z_len+1,:] .- c[2:kappa_z_len,:]) ./dz
    diff = (f_down .- f_up) ./ dz
    return diff
end
