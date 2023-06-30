
function rk4(ntemp, ptemp, ztemp, btemp, dtemp, otemp, prms, t)

    n, p, z, b, d, o = copy(ntemp), copy(ptemp), copy(ztemp), copy(btemp), copy(dtemp), copy(otemp)

    dNdt1, dPdt1, dZdt1, dBdt1, dDdt1, dOdt1 = model_functions(n, p, z, b, d, o, prms, t)
    
    track_n1 = n .+ prms.dt/2 .* dNdt1
    track_p1 = p .+ prms.dt/2 .* dPdt1
    track_z1 = z .+ prms.dt/2 .* dZdt1
    track_b1 = b .+ prms.dt/2 .* dBdt1
    track_d1 = d .+ prms.dt/2 .* dDdt1
    track_o1 = o .+ prms.dt/2 .* dOdt1

    N1tot = sum(track_p1) + sum(track_b1) + sum(track_z1) + sum(track_n1) + sum(track_d1)
    nan_or_inf(N1tot) && @error("Nan or Inf at t=$t: \n N1tot: $N1tot")
    
    dNdt2, dPdt2, dZdt2, dBdt2, dDdt2, dOdt2 = model_functions(track_n1, track_p1, track_z1, track_b1, track_d1, track_o1, prms, t)

    track_n2 = n .+ prms.dt/2 .* dNdt2
    track_p2 = p .+ prms.dt/2 .* dPdt2
    track_z2 = z .+ prms.dt/2 .* dZdt2
    track_b2 = b .+ prms.dt/2 .* dBdt2
    track_d2 = d .+ prms.dt/2 .* dDdt2
    track_o2 = o .+ prms.dt/2 .* dOdt2

    N2tot = sum(track_p2) + sum(track_b2) + sum(track_z2) + sum(track_n2) + sum(track_d2)
    nan_or_inf(N2tot) && @error("Nan or Inf at t=$t: \n N2tot: $N2tot") 
    
    dNdt3, dPdt3, dZdt3, dBdt3, dDdt3, dOdt3  = model_functions(track_n2, track_p2, track_z2, track_b2, track_d2, track_o2, prms, t)

    track_n3 = n .+ prms.dt .* dNdt3
    track_p3 = p .+ prms.dt .* dPdt3
    track_z3 = z .+ prms.dt .* dZdt3
    track_b3 = b .+ prms.dt .* dBdt3
    track_d3 = d .+ prms.dt .* dDdt3
    track_o3 = o .+ prms.dt .* dOdt3

    N3tot = sum(track_p3) + sum(track_b3) + sum(track_z3) + sum(track_n3) + sum(track_d3)
    nan_or_inf(N3tot) && @error("Nan or Inf at t=$t: \n N3tot: $N3tot") 
    
    dNdt4, dPdt4, dZdt4, dBdt4, dDdt4, dOdt4 = model_functions(track_n3, track_p3, track_z3, track_b3, track_d3, track_o3, prms, t)

    n .+= (dNdt1 .+ 2 .* dNdt2 .+ 2 .* dNdt3 .+ dNdt4) .* (prms.dt / 6)
    p .+= (dPdt1 .+ 2 .* dPdt2 .+ 2 .* dPdt3 .+ dPdt4) .* (prms.dt / 6)
    z .+= (dZdt1 .+ 2 .* dZdt2 .+ 2 .* dZdt3 .+ dZdt4) .* (prms.dt / 6)
    b .+= (dBdt1 .+ 2 .* dBdt2 .+ 2 .* dBdt3 .+ dBdt4) .* (prms.dt / 6)
    d .+= (dDdt1 .+ 2 .* dDdt2 .+ 2 .* dDdt3 .+ dDdt4) .* (prms.dt / 6)
    o .+= (dOdt1 .+ 2 .* dOdt2 .+ 2 .* dOdt3 .+ dOdt4) .* (prms.dt / 6)
     
    Ntot = sum(p) + sum(b) + sum(n) + sum(d) + sum(z)
    nan_or_inf(Ntot) && @error("Nan or Inf at t=$t: \n Ntot: $Ntot") 

    return n, p, z, b, d, o

end