# using Printf
# using Dates
# using NCDatasets
# using SparseArrays, LinearAlgebra


#make sum drop dimension automatically and set subnormals to zero to avoid slowdown
sumd(x,dims) = dropdims(sum(x,dims=dims),dims=dims)
set_zero_subnormals(true)


#TODO add a sinking (dsink) to remove some detritus or it wont equilibriate 
function run_NPZBD(prms, season)

    trec = prms.nt÷prms.nrec # frequency of recording
    start_time = now()

    # Generate empty arrays
    nrec1 = Int(prms.nrec + 1)
    track_n = Array{Float64, 3}(undef, prms.ngrid, prms.nn, nrec1) 
    track_p = Array{Float64, 3}(undef, prms.ngrid, prms.np, nrec1) 
    track_z = Array{Float64, 3}(undef, prms.ngrid, prms.nz, nrec1) 
    track_b = Array{Float64, 3}(undef, prms.ngrid, prms.nb, nrec1) 
    track_d = Array{Float64, 3}(undef, prms.ngrid, prms.nd, nrec1)
    track_o = Array{Float64, 3}(undef, prms.ngrid, 1, nrec1) 
    track_time = Array{Float64,1}(undef, nrec1)   
 
    track_n[:,:,1] .= prms.nIC
    track_p[:,:,1] .= prms.pIC
    track_z[:,:,1] .= prms.zIC
    track_b[:,:,1] .= prms.bIC
    track_d[:,:,1] .= prms.dIC
    track_o[:,:,1] .= prms.oIC
    track_time[1] = 0

    #--------------------------------------
    #Initial conditions at time t = 0
    ntemp = copy(prms.nIC) 
    ptemp = copy(prms.pIC) 
    ztemp = copy(prms.zIC) 
    btemp = copy(prms.bIC)
    dtemp = copy(prms.dIC) 
    otemp = copy(prms.oIC) 


    # @time for t = 1:prms.nt 
    for t = 1:prms.nt
        # Runge-Kutta 4th order 
        ntemp, ptemp, ztemp, btemp, dtemp, otemp = rk4(ntemp, ptemp, ztemp, btemp, dtemp, otemp, prms, t)

        if mod(t, trec) == 0

            track_n, track_p, track_z, track_b, track_d, track_o = update_tracking_arrs(track_n, track_p, track_z, track_b, track_d, track_o, track_time, 
                                                                                        ntemp, ptemp, ztemp, btemp, dtemp, otemp, t, trec, prms)
            println("Total N: ", sum(ptemp) + sum(btemp) + sum(ntemp) + sum(dtemp) + sum(ztemp))
            println("n = ", sum(ntemp))
            println("d = ", sum(dtemp))
            println("o = ", sum(otemp))
            
        end 

        # Nutrient pulsing routine 
        if prms.pulse != 1
            if season == 1
                if t % 1000 == 0 
                    println("PULSED at t=$t")
                    ntemp, dtemp = pulse_nutrients(ntemp, dtemp, prms, prms.pulse)      
                end
            else
                if t % 3000 == 0 
                    println("PULSED at t=$t")
                    ntemp, dtemp = pulse_nutrients(ntemp, dtemp, prms, prms.pulse)  
                end
            end
        end 

        #calculate bacteria uptake for last timepoint
        if t == prms.nt
            II, JJ = get_nonzero_axes(prms.CM)
            v = zeros(prms.nd, prms.nb) 
            uptake_b = zeros(prms.nd, prms.nb) 
        
            @inbounds for n = axes(II, 1)
                v[II[n],JJ[n]] = prms.umax_ij[II[n],JJ[n]] * prms.pen[JJ[n]] .* dtemp[II[n]] ./ (dtemp[II[n]] .+ prms.Km_ij[II[n],JJ[n]]) 
                uptake_b[II[n],JJ[n]] = v[II[n],JJ[n]] .* btemp[JJ[n]]
            end  
            
            end_time = now() 
            save_full_run(track_p, track_b, track_z, track_n, track_d, track_o, track_time, uptake_b, start_time, end_time, prms, season)
            save_endpoints(track_n, track_p, track_z, track_b, track_d, track_o, uptake_b, prms, season)

        end
    end 


    return ntemp, ptemp, ztemp, btemp, dtemp, otemp, track_time

end 


function model_functions(N, P, Z, B, D, O, prms, t)
    #TODO - figure out if there's any memory issues, is it increasing memory usage over time? seeems to be, figure out why
    #Transport
    dNdt = diffusion(N, prms.kappa_z, prms.dz)
    dPdt = diffusion(P, prms.kappa_z, prms.dz)
    dZdt = diffusion(Z, prms.kappa_z, prms.dz)
    dBdt = diffusion(B, prms.kappa_z, prms.dz)
    dDdt = diffusion(D, prms.kappa_z, prms.dz) .- advection(D, prms.wd, prms.dz)
    dOdt = diffusion(O, prms.kappa_z, prms.dz)

    d_gain_total = zeros(prms.ngrid)

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

    II, JJ = get_nonzero_axes(prms.CMp)

    for j = axes(II, 1)
        uptake = P[:,JJ[j]] .* prms.temp_fun .* prms.vmax_ij[II[j],JJ[j]] .* min.(N./ (N .+ prms.Kp_ij[II[j],JJ[j]]), prms.light ./ (prms.light .+ prms.K_I))
        dNdt += -uptake
        dOdt += uptake * prms.e_o
        dPdt[:,JJ[j]] += uptake 
    end

    return dPdt, dNdt, dOdt

end 


function bacteria_uptake(prms, B, D, dDdt, dBdt, dNdt, dOdt, t)

    II, JJ = get_nonzero_axes(prms.CM)

    for j = axes(II, 1)
        uptake = B[:,JJ[j]] .* prms.temp_fun .* prms.umax_ij[II[j],JJ[j]] .* D[:,II[j]] ./ (D[:,II[j]] .+ prms.Km_ij[II[j],JJ[j]])
        yield = prms.y_ij[II[j],JJ[j]]
        dDdt[:,II[j]] += -uptake
        dBdt[:,JJ[j]] += uptake .* yield
        dNdt += uptake .* (1 - yield)
        dOdt +=  -uptake .* (yield ./ prms.yo_ij[II[j],JJ[j]])
    end

    return dDdt, dBdt, dNdt, dOdt

end


function grazing(prms, P, B, Z, dZdt, dNdt, dPdt, dBdt, t)

    GrM = copy(prms.GrM)
    for k = 1:prms.nz
        if sum(GrM[k,1:prms.np]) > 0 
            dZdt, dNdt, dPdt = phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dPdt, k, t)
        end
        if sum(GrM[k,prms.np+1:end]) > 0 
            dZdt, dNdt, dBdt = bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dBdt, k, t)
        end
    end

    return dZdt, dNdt, dPdt, dBdt

end

        function phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dPdt, k, t)

            prey = sum(GrM[k,1:prms.np]' .*P, dims=2)
            gp = prms.g_max[k] .* prey ./ (prey .+ prms.K_g[k])
            dZdt[:,k] += prms.γ[k] .* gp .* Z[:,k]
            dNdt += (1 - prms.γ[k]) .* gp .* Z[:,k]
            dPdt += -gp .* Z[:,k] .* GrM[k,1:prms.np]' .* P ./ prey 
            
            return dZdt, dNdt, dPdt

        end

        function bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dBdt, k, t)

            prey = sum(GrM[k,prms.np+1:end]' .*B, dims=2)
            gb = prms.g_max[k] .* prey ./ (prey .+ prms.K_g[k])
            dZdt[:,k] += prms.γ[k] .* gb .* Z[:,k]
            dNdt += (1 - prms.γ[k]) .* gb .* Z[:,k]
            dBdt +=  -gb .* Z[:,k] .* GrM[k,prms.np+1:end]' .* B ./ prey

            return dZdt, dNdt, dBdt

        end


function phyto_mortality(prms, P, dPdt, d_gain_total, t)

    pmort = (transpose(prms.m_lp) .+ transpose(prms.m_qp) .* P) .* P
    dPdt += -pmort
    d_gain_total += sum(pmort, dims=2)

    return dPdt, d_gain_total

end


function bacterial_mortality(prms, B, dBdt, d_gain_total, t)

    bmort = (transpose(prms.m_lb) .+ transpose(prms.m_qb) .* B) .* B
    dBdt += -bmort
    d_gain_total += sum(bmort, dims=2)

    return dBdt, d_gain_total

end


function zoo_mortality(prms, Z, dZdt, d_gain_total, t)

    zmort = (transpose(prms.m_lz) .+ transpose(prms.m_qz) .* Z) .* Z
    dZdt += -zmort
    d_gain_total += sum(zmort, dims=2)

    return dZdt, d_gain_total

end

function change_in_o2(prms, O, dOdt, t)

    dOdt[1:Int(prms.ml_boxes)] += prms.koverh .* (prms.o2_sat .- O[1:Int(prms.ml_boxes)])
    dOdt += prms.t_o2relax .* (prms.o2_deep .- O) #relaxation at depth (lateral flux)

    return dOdt

end

function distribute_d_gain(prms, dDdt, d_gain_total, t)

    dDdt += d_gain_total .* transpose(prms.prob_generate_d)

    return dDdt

end


function get_nonzero_axes(M)

    Cs = sparse(M)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 


