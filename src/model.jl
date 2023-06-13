using Printf
using Dates
using NCDatasets
using SparseArrays, LinearAlgebra


#make sum drop dimension automatically and set subnormals to zero to avoid slowdown
sumd(x,dims) = dropdims(sum(x,dims=dims),dims=dims)
set_zero_subnormals(true)


#TODO add a sinking (dsink) to remove some detritus or it wont equilibriate 
function run_NPZBD(prms)

    start_time = now()

    # dt = prms.dt
    nt = Int(tt/dt)
    trec = nt÷nrec # frequency of recording
    fsaven = print_info(start_time, prms, nt)

    # Generate empty arrays
    nrec1 = Int(nrec + 1)
    track_n = Array{Float64, 3}(undef, ngrid, nn, nrec1) 
    track_p = Array{Float64, 3}(undef, ngrid, np, nrec1) 
    track_z = Array{Float64, 3}(undef, ngrid, nz, nrec1) 
    track_b = Array{Float64, 3}(undef, ngrid, nb, nrec1) 
    (track_d, SW_all) = (Array{Float64, 3}(undef, ngrid, nd, nrec1) for _ in 1:2)
    track_o = Array{Float64, 3}(undef, ngrid, 1, nrec1) 
    track_time = Array{Float64,1}(undef, nrec1)   
 
    track_n[:,:,1] .= nIC
    track_p[:,:,1] .= pIC
    track_z[:,:,1] .= zIC
    track_b[:,:,1] .= bIC
    track_d[:,:,1] .= dIC
    track_o[:,:,1] .= oIC
    SW_all[:,:,1] .= NaN
    track_time[1] = 0

    #--------------------------------------
    #Initial conditions at time t = 0
    ntemp = copy(nIC) 
    ptemp = copy(pIC) 
    ztemp = copy(zIC) 
    btemp = copy(bIC)
    dtemp = copy(dIC) 
    otemp = copy(oIC) 


    for t = 1:nt 

        # Runge-Kutta 4th order 
        ntemp, ptemp, ztemp, btemp, dtemp, otemp = rk4(ntemp, ptemp, ztemp, btemp, dtemp, otemp, prms, t)

        if mod(t, trec)==0

            track_n, track_p, track_z, track_b, track_d, track_o = update_tracking_arrs(track_n, track_p, track_z, track_b, track_d, track_o, track_time, 
                                                                                        ntemp, ptemp, ztemp, btemp, dtemp, otemp, t, trec, prms)
            println("Total N: ", sum(ptemp) + sum(btemp) + sum(ntemp) + sum(dtemp) + sum(ztemp))
        
        end 

        # nutrient pulse every 7 days
        # if t % 700 == 0
        #     pulse = nutrient_pulse()
        #     ntemp = ntemp .+ pulse           
        # end

        #calculate uptake and v for last timepoint
        if t == nt

            II, JJ = get_nonzero_axes(prms)
            v = zeros(nd,nb) 
            uptake = zeros(nd,nb) 
        
            @inbounds for n = axes(II, 1)
                v[II[n],JJ[n]] = vmax_ij[II[n],JJ[n]] * pen[JJ[n]] .* dtemp[II[n]] ./ (dtemp[II[n]] .+ Km_ij[II[n],JJ[n]]) 
                uptake[II[n],JJ[n]] = v[II[n],JJ[n]] .* btemp[JJ[n]]
            end  
            
            end_time = now() 

            savetoNC(fsaven, track_p, track_b, track_z, track_n, track_d, track_o, track_time, v, uptake, start_time, end_time, params)

        end
    end 


    return ntemp, ptemp, ztemp, btemp, dtemp, otemp, track_time, fsaven

end 


function model_functions(N, P, Z, B, D, O, prms, t)

    #Transport
    dNdt = diffusion(N, kappa_z, dz)
    dPdt = diffusion(P, kappa_z, dz)
    dZdt = diffusion(Z, kappa_z, dz)
    dBdt = diffusion(B, kappa_z, dz)
    dDdt = diffusion(D, kappa_z, dz) .- advection(D, wd, dz)
    dOdt = diffusion(O, kappa_z, dz)

    d_gain_total = zeros(ngrid)

    # phyto uptake 
    #TODO phyto assimilation efficiency and contribution to dDdt?
    dPdt, dNdt, dOdt = phyto_uptake(prms, N, P, dNdt, dPdt, dOdt, t)

    # bacteria uptake
    dDdt, dBdt, dNdt, dOdt = bacteria_uptake(prms, B, D, dDdt, dBdt, dNdt, dOdt, t)

    # grazing
    dZdt, dNdt, dPdt, dBdt = grazing(prms, P, B, Z, dZdt, dNdt, dPdt, dBdt, t)

    #phytoplankton mortality 
    dPdt, d_gain_total = phyto_mortality(prms, P, dPdt, d_gain_total, t)

    #bacterial mortality 
    dBdt, d_gain_total = bacterial_mortality(prms, B, dBdt, d_gain_total, t)
    
    #zooplankton mortality 
    dZdt, d_gain_total = zoo_mortality(prms, Z, dZdt, d_gain_total, t)

    #split accumulated OM into nd pools according to probability of generation
    dDdt = distribute_d_gain(prms, dDdt, d_gain_total, t)

    #oxygen
    dOdt = change_in_o2(prms, O, dOdt, t)


    return dNdt, dPdt, dZdt, dBdt, dDdt, dOdt

end 


function phyto_uptake(prms, N, P, dNdt, dPdt, dOdt, t)
    #NOTE > only works for nn = 1 at present. Implement CM for p_i on n_i and adapt loops accordingly
    
    for i = 1:np
        uptake = P[:,i] .* temp_fun .* umax_p[i] .* min.(N ./ (N .+ K_n[i]), light ./ (light .+ K_I[i]))
        dNdt += -uptake
        dOdt += uptake * e_o
        dPdt[:,i] += uptake 
    end

    return dPdt, dNdt, dOdt

end 


function bacteria_uptake(prms, B, D, dDdt, dBdt, dNdt, dOdt, t)

    II, JJ = get_nonzero_axes(prms)

    for j = axes(II, 1)
        uptake = vmax_ij[II[j],JJ[j]] * pen[JJ[j]] * D[II[j]] / (D[II[j]] + Km_ij[II[j],JJ[j]]) * B[JJ[j]]
        uptake = B[:,JJ[j]] .* temp_fun .* vmax_ij[II[j],JJ[j]] .* D[:,II[j]] ./ (D[:,II[j]] .+ Km_ij[II[j],JJ[j]]) 
        dDdt[:,II[j]] += -uptake
        dBdt[:,JJ[j]] += uptake .* y_ij[II[j],JJ[j]]
        dNdt += uptake .* (1 - y_ij[II[j],JJ[j]])
        dOdt +=  -uptake .* y_ij[II[j],JJ[j]] ./ yo_ij[II[j],JJ[j]]
    end

    return dDdt, dBdt, dNdt, dOdt

end


function grazing(prms, P, B, Z, dZdt, dNdt, dPdt, dBdt, t)

    GrM = copy(prms.GrM)
    for k = 1:prms.nz
        if sum(GrM[k,1:np+1]) > 0 
            dZdt, dNdt, dPdt = phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dPdt, k, t)
        end
        if sum(GrM[k,np+1:end]) > 0 
            dZdt, dNdt, dBdt = bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dBdt, k, t)
        end
    end

    return dZdt, dNdt, dPdt, dBdt

end

        function phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dPdt, k, t)

            prey = sum(GrM[k,1:np+1]' .*P, dims=2)
            gp = g_max[k] .* prey ./ (prey .+ K_g[k])
            dZdt[:,k] += γ[k] .* gp .* Z[:,k]
            dNdt += (1 - γ[k]) .* gp .* Z[:,k]
            dPdt += -gp .* Z[:,k] .* GrM[k,1:np]' .* P ./ prey 
            
            return dZdt, dNdt, dPdt

        end

        function bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dBdt, k, t)

            prey = sum(GrM[k,np+1:end]' .*B, dims=2)
            gb = g_max[k] .* prey ./ (prey .+ K_g[k])
            dZdt[:,k] += γ[k] .* gb .* Z[:,k]
            dNdt += (1 - γ[k]) .* gb .* Z[:,k]
            dBdt +=  -gb .* Z[:,k] .* GrM[k,np+1:end]' .* B ./ prey

            return dZdt, dNdt, dBdt

        end


function phyto_mortality(prms, P, dPdt, d_gain_total, t)

    pmort = (transpose(m_lp) .+ transpose(m_qp) .* P) .* P
    dPdt += -pmort
    d_gain_total += sum(pmort, dims=2)

    return dPdt, d_gain_total

end


function bacterial_mortality(prms, B, dBdt, d_gain_total, t)

    bmort = (transpose(m_lb) .+ transpose(m_qb) .* B) .* B
    dBdt += -bmort
    d_gain_total += sum(bmort, dims=2)

    return dBdt, d_gain_total

end


function zoo_mortality(prms, Z, dZdt, d_gain_total, t)

    zmort = (transpose(m_lz) .+ transpose(m_qz) .* Z) .* Z
    dZdt += -zmort
    d_gain_total += sum(zmort, dims=2)

    return dZdt, d_gain_total

end

function change_in_o2(prms, O, dOdt, t)

    dOdt[1:Int(ml_boxes)] += koverh .* (o2_sat .- O[1:Int(ml_boxes)])
    dOdt += t_o2relax .* (o2_deep .- O) #relaxation at depth (lateral flux)

    return dOdt

end

function distribute_d_gain(prms, dDdt, d_gain_total, t)

    dDdt += d_gain_total .* transpose(prob_generate_d)

    return dDdt

end


function get_nonzero_axes(prms)

    Cs = sparse(prms.CM)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 


function nutrient_pulse()

    pulse = 5.0

    return pulse

end
