using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames, CSV
using SparseArrays, LinearAlgebra


#--------------------------------------------------------------------------------------
# R* Bacteria on Detritus
#--------------------------------------------------------------------------------------
function get_rstar_B(B, Z, ds, nb, nz, season)
    
    mort_b = b_mortality(B, ds, nb)
    grz_b = b_grazing(B, Z, ds, nb, nz)
    loss_b = b_loss(mort_b, grz_b, nb)
    RstarB_ij = Rstar(loss_b, ds, season)

    return RstarB_ij

end


function b_mortality(B, ds, nb)

    mort_b = Any[]
    for i in range(1, nb)
        push!(mort_b, (ds["m_lb"][i] .+ ds["m_qb"][i] .* B[:,i]))
    end

    return mort_b

end


function b_grazing(B, Z, ds, nb, nz)
    #------- for 1N 4P 3Z 7B 4D
    if nb==7
        GrM = ds["GrM"][:]
        grazing = Any[]
        g_max = 1.0
        K_g = ds["K_g"][3]
        K_g_POM = ds["K_g"][2]

        prey_POM = GrM[2,5]' .*B[:,1]
        gb_POM = g_max .* prey_POM ./ (prey_POM .+ K_g_POM)
        # grz_POM = gb_POM .* Z[:,2] .* GrM[2,5]'  .* B[:,1] ./ prey_POM
        grz_POM = (gb_POM .* Z[:,2] .* GrM[2,5]'  .* B[:,1] ./ prey_POM) ./ B[:,1] 
        push!(grazing, grz_POM)

        prey = GrM[3,6:end]' .*B[:,2:end]
        for i in range(1, 6)
            gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
            # grz_i = gb_i .* Z[:,3] .* GrM[3,i+5]'  .* B[:,i+1] ./ prey[:,i]
            grz_i = (gb_i .* Z[:,3] .* GrM[3,i+5]' ./ prey[:,i]) 
            push!(grazing, grz_i)
        end
    
        return grazing
    
    else

        #----- for 1N 2P 2Z 2B 2D
        if nz == 2
            GrM = ds["GrM"][:]
            grazing = Any[]
            g_max = 1.0
            K_g = ds["K_g"][2]

            prey = GrM[2,3:end]' .*B[:,1:end]
            for i in range(1, nb)
                gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
                grz_i = (gb_i .* Z[:,2] .* GrM[2,i+2]' ./ prey[:,i])
                push!(grazing, grz_i)
            end

        #----- for 1N 1P 3Z 2B 2D
        elseif nz == 3
            GrM = ds["GrM"][:]
            grazing = Any[]
            g_max = 1.0
            K_g = [ds["K_g"][2], ds["K_g"][3]]

            prey = 1.0 .*B[:,1:2]
            for i in range(1, nb)
                gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g[i])
                grz_i = (gb_i .* Z[:,i+1] .* 1.0 ./ prey[:,i]) 
                push!(grazing, grz_i)
            end
        else
        end
    end

    return grazing

end


function b_loss(mortality, grazing, nb)

    loss = Any[]
    for i in range(1, nb)
        push!(loss, mortality[i] .+ grazing[i])
    end

    return loss
end

function get_temp_mod(season)
    #fit to SPOT data (approx 20 to 4, approx 16 to 4)
    if season == "Win"
        temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/win_temp_mod.csv", DataFrame)
    else
        temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/sum_temp_mod.csv", DataFrame)
    end

    return Matrix(temp_mod)
end

function Rstar(loss, ds, season)

    vmax_ij = ds["vmax_ij"][:]
    Km_ij = ds["Km_ij"][:]
    yield = ds["y_ij"][:]
    temp_mod = get_temp_mod(season)

    II, JJ = get_nonzero_axes(ds["CM"][:])
    RS = Any[]
    for j = axes(II, 1)
        push!(RS, Km_ij[II[j],JJ[j]] .* loss[j] ./ (yield[II[j],JJ[j]] .* vmax_ij[II[j],JJ[j]] .* temp_mod .- loss[j]))
    end

    for i in range(1, length(RS))
        RS[i] = check_for_negatives(RS[i])
    end

    return RS
end


#--------------------------------------------------------------------------------------
# R* Phyto on Nutrients
#--------------------------------------------------------------------------------------

function get_rstar_P(P, Z, ds, n, season)
    
    mort_p = p_mortality(P, ds, n)
    grz_p = p_grazing(P, Z, ds, n)
    loss_p = p_loss(mort_p, grz_p, n)
    RstarP_i = RstarP(loss_p, ds, n, season)

    return RstarP_i

end

function p_mortality(P, ds, n)

    m_lp = 1e-1
    m_qp = 1e-1
    mort_p = Any[]
    for i in range(1, n)
        push!(mort_p, (m_lp .+ m_qp .* P[:,i]))
    end

    return mort_p

end

function p_grazing(P, Z, ds, n)

    GrM = ds["GrM"][:]
    grazing = Any[]
    g_max = 1.0
    K_g = ds["K_g"][1]

    prey = GrM[1,1:n]' .*P[:,:]
    for i in range(1, n)
        gp_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
        grz_i = (gp_i .* Z[:,1] .* GrM[1,i]') ./ prey[:,i]
        push!(grazing, grz_i)
    end

    return grazing

end

function p_loss(mortality, grazing, n)

    loss = Any[]
    for i in range(1, n)
        push!(loss, mortality[i] .+ grazing[i])
    end

    return loss
end

function RstarP(loss, ds, n, season)

    umax_ij = ds["umax_ij"][:]
    Kp_ij = ds["Kp_ij"][:]
    temp_mod = get_temp_mod(season)

    RS = Any[]
    for i =  range(1, n)
        push!(RS, Kp_ij[i] .* loss[i] ./ (umax_ij[i] .* temp_mod .- loss[i]))
    end

    for i in range(1, length(RS))
        RS[i] = check_for_negatives(RS[i])
    end

    return RS
end

#--------------------------------------------------------------------------------------
# R* Z on B & P
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
# Plotting
#--------------------------------------------------------------------------------------
function plot_rstar_7B4D(rstar_w, rstar_s, Dw, Ds, Bw, Bs, ds)
    # B1 eats D1 (POM)  |  B2 & B7 eat D2  |  B3 & B6 eat D3  |  B4 & B5 eat D4

    H = 500
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]
    lg=:bottomright
    tfs=18

    #--------------------------------------------------------------------------------------------
    # Winter
    #--------------------------------------------------------------------------------------------
    p1 = plot(rstar_w[1][1:50], -zc, lw=4, lc="red3", label=L" B1", ylabel="Depth (m)", 
    xrotation=45, xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10)
    plot!(Dw[1:50, 1], -zc, lw=3, lc="blue3", linestyle=:dot,label=L" D1", alpha=0.4)

    p2 = plot(rstar_w[3][1:50], -zc, lw=4, lc="darkorange", label="", 
    xrotation=45, xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    plot!(rstar_w[6][1:50], -zc, lw=4, lc="magenta2",  label="")
    plot!(Dw[1:50, 3], -zc, lw=4, lc="black", linestyle=:dot, alpha=0.6, label=L" D3")

    p3 = plot(Bw[1:50, 3], -zc, lw=4, lc="darkorange", label=L" B3", 
    xrotation=45, xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), xlims=(-0.01, 0.2))
    plot!(Bw[1:50, 6], -zc, lw=4, lc="magenta2",  label=L" B6", alpha=0.8)

    p4 = plot(rstar_w[4][1:50], -zc, lw=4, lc="gold4", label="", 
    xrotation=45, xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    plot!(rstar_w[5][1:50], -zc, lw=4, lc="lightcyan4", label="")
    plot!(Dw[1:50, 4], -zc, lw=4, lc="purple", linestyle=:dot, alpha=0.6, label=L" D4")

    p5 = plot(Bw[1:50, 4], -zc, lw=4, lc="gold4", label=L" B4", 
    xrotation=45, xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), xlims=(-0.01, 0.2))
    plot!(Bw[1:50, 5], -zc, lw=4, lc="lightcyan4", label=L" B5", alpha=0.8)

    #--------------------------------------------------------------------------------------------
    # Summer
    #--------------------------------------------------------------------------------------------
    p6 = plot(rstar_s[1][1:50], -zc, lw=4, lc="red3", label=L" B1", ylabel="Depth (m)", 
    xrotation=45, xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10)
    plot!(Ds[1:50, 1], -zc, lw=3, lc="blue3", linestyle=:dot,label=L" D1", alpha=0.4)

    p7 = plot(rstar_s[3][1:50], -zc, lw=4, lc="darkorange", label="", 
    xrotation=45, xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    plot!(rstar_s[6][1:50], -zc, lw=4, lc="magenta2", label="")
    plot!(Ds[1:50, 3], -zc, lw=4, lc="black", linestyle=:dot, alpha=0.6, label=L" D3")

    p8 = plot(Bs[1:50, 3], -zc, lw=4, lc="darkorange", label=L" B3", 
    xrotation=45, xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), xlims=(-0.01, 0.2))
    plot!(Bs[1:50, 6], -zc, lw=4, lc="magenta2", label=L" B6", alpha=0.8)

    p9 = plot(rstar_s[4][1:50], -zc, lw=4, lc="gold4", label="", 
    xrotation=45, xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    plot!(rstar_s[5][1:50], -zc, lw=4, lc="lightcyan4", label="")
    plot!(Ds[1:50, 4], -zc, lw=4, lc="purple", linestyle=:dot, alpha=0.6, label=L" D4")

    p10 = plot(Bs[1:50, 4], -zc, lw=4, lc="gold4", label=L" B4", 
    xrotation=45, xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), xlims=(-0.01, 0.2))
    plot!(Bs[1:50, 5], -zc, lw=4, lc="lightcyan4", label=L" B5", alpha=0.8)
    
    
    combined = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
        fg_legend = :transparent,
        layout = (2,5),
        size=(900,700),
    )

    savefig(combined, "/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rstar_7B4D.png")

end


function plot_rstar_2B2D(rstar_w, rstar_s, Dw, Ds, Bw, Bs, ds, np)
    # B1 eats D1 (POM)  |  B2 eats D2 (DOM)  
    # zc = get_zc(ds, 500)
    H = 500
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]
    lg=:bottomright
    tfs=12

    #--------------------------------------------------------------------------------------------
    # Winter
    #--------------------------------------------------------------------------------------------
    p1 = plot(rstar_w[1][1:50], -zc, lw=4, lc="red3", ylabel="Depth (m)", label="", xrotation=45, 
    xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, alpha=0.7)
    plot!(Dw[1:50, 1], -zc, lw=3, lc="blue3", linestyle=:dot, alpha=0.7, label=L" D1")

    p2 = plot(rstar_w[2][1:50], -zc, lw=4, lc="darkorange", label="",xrotation=45, 
    xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), alpha=0.7)
    plot!(Dw[1:50, 2], -zc, lw=3, lc="black", linestyle=:dot, alpha=0.7, label=" D2")

    p3 = plot(Bw[1:50, 1], -zc, lw=4, lc="red3", ylabel="", label=" B1",xrotation=45, 
    xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), alpha=0.7)
    plot!(Bw[1:50, 2], -zc, lw=4, lc="darkorange", ylabel="", label=" B2", xrotation=45, 
    xguidefontsize=12, border=:box, legend=lg, alpha=0.7)

    #--------------------------------------------------------------------------------------------
    # Summer
    #--------------------------------------------------------------------------------------------
    p4 = plot(rstar_s[1][1:50], -zc, lw=4, lc="red3", ylabel="Depth (m)", label="", xrotation=45, 
    xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, alpha=0.7)
    plot!(Ds[1:50, 1], -zc, lw=3, lc="blue3", linestyle=:dot, alpha=0.7, label=L" D1")

    p5 = plot(rstar_s[2][1:50], -zc, lw=4, lc="darkorange", label="",xrotation=45, 
    xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), alpha=0.7)
    plot!(Ds[1:50, 2], -zc, lw=3, lc="black", linestyle=:dot, alpha=0.7, label=" D2")

    p6 = plot(Bs[1:50, 1], -zc, lw=4, lc="red3", ylabel="", label=" B1",xrotation=45, 
    xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), alpha=0.7)
    plot!(Bs[1:50, 2], -zc, lw=4, lc="darkorange", ylabel="", label=" B2", xrotation=45, 
    xguidefontsize=12, border=:box, legend=lg, alpha=0.7)

    combined = plot(p1, p2, p3, p4, p5, p6,
        fg_legend = :transparent,
        fontfamily="Computer Modern",
        layout = (2,3),
        size=(600,500),
    )

    savefig(combined, "/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rstarB_$(np)P2B2D.png")

end



function plot_rstar_P(rstar_w, rstar_s, Pw, Ps, Nw, Ns, ds, np)

    H = 200
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]
    lg=:bottomright
    tfs=12

    if np == 1
        p1 = plot(rstar_w[1][1:20], -zc, lw=4, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(Nw[1:20, 1], -zc, lw=4, lc="blue3", linestyle=:dot, label=L" N", alpha=0.4, legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=4, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))

        p3 = plot(rstar_s[1][1:20], -zc, lw=4, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(Ns[1:20, 1], -zc, lw=4, lc="blue3", linestyle=:dot,label=L" N", alpha=0.4, legend=:bottom)

        p4 = plot(Ps[1:20, 1], -zc, lw=4, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))

    
    elseif np == 2
        p1 = plot(rstar_w[1][1:20], -zc, lw=4, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(rstar_w[2][1:20], -zc, lw=4, lc="darkgreen",  label="")
        plot!(Nw[1:20, 1], -zc, lw=4, lc="blue3", linestyle=:dot, label=L" N", alpha=0.4, legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=4, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))
        plot!(Pw[1:20, 2], -zc, lw=4, lc="darkgreen",  label=L" P2")

        p3 = plot(rstar_s[1][1:20], -zc, lw=4, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(rstar_s[2][1:20], -zc, lw=4, lc="darkgreen",  label="")
        plot!(Ns[1:20, 1], -zc, lw=4, lc="blue3", linestyle=:dot,label=L" N", alpha=0.4, legend=:bottom)

        p4 = plot(Ps[1:20, 1], -zc, lw=4, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))
        plot!(Ps[1:20, 2], -zc, lw=4, lc="darkgreen",  label=L" P2")

    elseif np == 4
        p1 = plot(rstar_w[1][1:20], -zc, lw=4, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(rstar_w[2][1:20], -zc, lw=4, lc="darkgreen",  label="")
        plot!(rstar_w[3][1:20], -zc, lw=4, lc="red4",  label="")
        plot!(rstar_w[4][1:20], -zc, lw=4, lc="cyan4",  label="")
        plot!(Nw[1:20, 1], -zc, lw=4, lc="blue3", linestyle=:dot, label=L" N", alpha=0.4, legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=4, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))
        plot!(Pw[1:20, 2], -zc, lw=4, lc="darkgreen",  label=L" P2")
        plot!(Pw[1:20, 3], -zc, lw=4, lc="red4",  label=L" P3")
        plot!(Pw[1:20, 4], -zc, lw=4, lc="cyan4",  label=L" P4")

        p3 = plot(rstar_s[1][1:20], -zc, lw=4, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(rstar_s[2][1:20], -zc, lw=4, lc="darkgreen",  label="")
        plot!(rstar_s[3][1:20], -zc, lw=4, lc="red4",  label="")
        plot!(rstar_s[4][1:20], -zc, lw=4, lc="cyan4",  label="")
        plot!(Ns[1:20, 1], -zc, lw=4, lc="blue3", linestyle=:dot,label=L" N", alpha=0.4, legend=:bottom)

        p4 = plot(Ps[1:20, 1], -zc, lw=4, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))
        plot!(Ps[1:20, 2], -zc, lw=4, lc="darkgreen",  label=L" P2")
        plot!(Ps[1:20, 3], -zc, lw=4, lc="red4",  label=L" P3")
        plot!(Ps[1:20, 4], -zc, lw=4, lc="cyan4",  label=L" P4")
        
    else
    end
    

    combined = plot(p1, p2, p3, p4,
    fg_legend = :transparent,
    layout = (2,2),
    size=(450,600),
    plot_title=L"R* P_N"
    )

    savefig(combined,"/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rstarP_1N$(np)P.png")

end


#--------------------------------------------------------------------------------------
# Utility Functions
#--------------------------------------------------------------------------------------

function get_endpoints(ds, vars)

    out = Vector{Any}()

    for v in vars
        append!(out, [ds["$v"][:,:,end]])
    end

    return out[1], out[2], out[3], out[4], out[5]

end

function get_nonzero_axes(M)

    Cs = sparse(M)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 

function cut_off(ds, n)

    dss = copy(ds)
    co = 10^-6
    for i in range(1, n)
        dss[:, i] .= ifelse.(dss[:, i] .< co, co, dss[:, i])
    end

    return dss

end

function extinct(ds, n)

    dss = copy(ds)
    ex = 10^-6
    for i in range(1, n)
        dss[:, i] .= ifelse.(dss[:, i] .== ex, 0.0, dss[:, i])
    end

    return dss

end

function check_for_negatives(RS)

    for i in eachindex(RS)
        RS[i] = ifelse(RS[i] < 0, NaN, RS[i])
    end

    return RS

end

function get_zc(ds, H=890)

    dz = 10
    zc = [dz/2:dz:(H-dz/2)]

    return zc
end


#--------------------------------------------------------------------------------------
# ---------- RUN
#--------------------------------------------------------------------------------------


# -----------------------------------------------------------------------
# 1N 4P 3Z 7B 4D
# -----------------------------------------------------------------------
winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230827_1345.nc")
summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230827_1710.nc")
CMd = display(winter["CM"][:])
GrMd = display(winter["GrM"][:])

# Get endpoints and set all biomass values below cut-off (10^-6) equal to cut-off
Nw, Pw, Zw, Bw, Dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
Ns, Ps, Zs, Bs, Ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

Bww = cut_off(Bw, 7)
Bss = cut_off(Bs, 7)
Pww = cut_off(Pw, 4)
Pss = cut_off(Ps, 4)

# Calculate Rstar B
rstarB_ijw = get_rstar_B(Bww, Zw, winter, 7, 3, "Win")
rstarB_ijs = get_rstar_B(Bss, Zs, summer, 7, 3, "Sum")
rstarP_w = get_rstar_P(Pww, Zw, winter, 4, "Win")
rstarP_s = get_rstar_P(Pss, Zs, summer, 4, "Sum")

Bwe = extinct(Bww, 7)
Bse = extinct(Bss, 7)
Pwe = extinct(Pww, 4)
Pse = extinct(Pss, 4)

plot_rstar_7B4D(rstarB_ijw, rstarB_ijs, Dw, Ds, Bwe, Bse, winter)
plot_rstar_P(rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, winter,  4)
println("done")


# -----------------------------------------------------------------------
# 1N 2P 2Z 2B 2D
# -----------------------------------------------------------------------
winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230905_2005.nc")
summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230905_2301.nc")
CMd = display(winter["CM"][:])
GrMd = display(winter["GrM"][:])

# Get endpoints 
Nw, Pw, Zw, Bw, Dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
Ns, Ps, Zs, Bs, Ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

# set all biomass values below cut-off (10^-6) equal to cut-off
Bww = cut_off(Bw, 2)
Bss = cut_off(Bs, 2)
Pww = cut_off(Pw, 2)
Pss = cut_off(Ps, 2)

# Calculate Rstar 
rstarB_ijw = get_rstar_B(Bww, Zw, winter, 2, 2, "Win")
rstarB_ijs = get_rstar_B(Bss, Zs, summer, 2, 2, "Sum")
rstarP_w = get_rstar_P(Pww, Zw, winter, 2, "Win")
rstarP_s = get_rstar_P(Pss, Zs, summer, 2, "Sum")

# set low biomass to zero
Bwe = extinct(Bww, 2)
Bse = extinct(Bss, 2)
Pwe = extinct(Pww, 2)
Pse = extinct(Pss, 2)


# Plot
plot_rstar_2B2D(rstarB_ijw, rstarB_ijs, Dw, Ds, Bwe, Bse, winter, 2)
plot_rstar_P(rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, winter,  2)
println("done")


# -----------------------------------------------------------------------
# 1N 1P 3Z 2B 2D
# -----------------------------------------------------------------------
winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230906_1351.nc")
summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230906_1413.nc")
CMd = display(winter["CM"][:])
GrMd = display(winter["GrM"][:])

# Get endpoints
Nw, Pw, Zw, Bw, Dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
Ns, Ps, Zs, Bs, Ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

# set all biomass values below cut-off (10^-6) equal to cut-off
Bww = cut_off(Bw, 2)
Bss = cut_off(Bs, 2)
Pww = cut_off(Pw, 1)
Pss = cut_off(Ps, 1)

# Calculate Rstar 
rstarB_ijw = get_rstar_B(Bww, Zw, winter, 2, 3, "Win")
rstarB_ijs = get_rstar_B(Bss, Zs, summer, 2, 3, "Sum")
rstarPw = get_rstar_P(Pww, Zw, winter, 1, "Win")
rstarPs = get_rstar_P(Pss, Zs, summer, 1, "Sum")

Bwe = extinct(Bww, 2)
Bse = extinct(Bss, 2)
Pwe = extinct(Pww, 1)
Pse = extinct(Pss, 1)

# Plot
plot_rstar_2B2D(rstarB_ijw, rstarB_ijs, Dw, Ds, Bwe, Bse, winter, 1)
plot_rstar_P(rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, winter,  1)
println("done")



# function get_rstar_B(B, Z, ds)
    
#     mort_b = b_mortality(B, ds)
#     grz_b = b_grazing(B, Z, ds)
#     loss_b = b_loss(mort_b, grz_b)
#     RstarB_ij = Rstar(loss_b, ds)

#     return RstarB_ij

# end

# function b_mortality(B, ds)

#     mort_b = Any[]
#     for i in range(1, 7)
#         push!(mort_b, (ds["m_lb"][i] .+ ds["m_qb"][i] .* B[:,i]))
#     end

#     return mort_b

# end


# function b_grazing(B, Z, ds)

#     GrM = ds["GrM"][:]
#     grazing = Any[]
#     g_max = 1.0
#     K_g = ds["K_g"][3]
#     K_g_POM = ds["K_g"][2]

#     prey_POM = GrM[2,5]' .*B[:,1]
#     gb_POM = g_max .* prey_POM ./ (prey_POM .+ K_g_POM)
#     grz_POM = gb_POM .* Z[:,2] .* GrM[2,5]'  .* B[:,1] ./ prey_POM
#     push!(grazing, grz_POM)

#     prey = GrM[3,6:end]' .*B[:,2:end]
#     for i in range(1, 6)
#         gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
#         grz_i = gb_i .* Z[:,3] .* GrM[3,i+5]' ./ prey[:,i]
#         push!(grazing, grz_i)
#     end

#     return grazing

# end

# function b_loss(mortality, grazing)

#     loss = Any[]
#     for i in range(1, 7)
#         push!(loss, mortality[i] .+ grazing[i])
#     end

#     return loss
# end

# function Rstar(loss, ds)

#     II, JJ = get_nonzero_axes(ds["CM"][:])
#     vmax_ij = ds["vmax_ij"][:]
#     Km_ij = ds["Km_ij"][:]

#     RS = Any[]
#     for j = axes(II, 1)
#         push!(RS, Km_ij[II[j],JJ[j]] .* loss[j] ./ (vmax_ij[II[j],JJ[j]] .- loss[j]))
#     end

#     return RS
# end