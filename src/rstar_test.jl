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
    RstarB_ij = Rstar(loss_b, ds, season, nb)

    return RstarB_ij

end


function b_mortality(B, ds, nb)

    mort_b = Any[]
    for i in range(1, nb)
        push!(mort_b, (ds["m_lb"][i] .+ ds["m_qb"][i] .* B[:,i]))
    end

    return mort_b

end

# function get_prey(GrM, B)
#     B_dom = B[:,2:end]
#     GrM_dom = GrM[4:end,10:end]

#     for (i, row) in enumerate(eachrow(GrM_dom))
#         for (j, col) in enumerate(eachcol(GrM_dom))
#             if GrM_dom[i, j] > 0
#                 prey = GrM[i, j] .* B[:,j]
#             end
#         end
#     end
# end

function b_grazing(B, Z, ds, nb, nz)
    #NOTE this all needs to be generalized / scaled !! 
    GrM = ds["GrM"][:]
    grazing = Any[]
    g_max = 1.0
    K_g = 1.0

    #------- for 1N 8P 6Z 13B 5D
    if nb == 13
        # NOTE length(B[:,1]) = prms.ngrid
        grazing_zi = zeros(Float64, length(B[:,1]), nb) 
        np = 8
        for z in range(1, nz)
            if sum(GrM[z,np+1:end]) > 0 
                prey = GrM[z,np+1:end]' .*B[:,1:end]
                g = g_max .* prey ./ (prey .+ K_g)
                grazing_zi += (g .* Z[:,z] .* GrM[z,np+1:end]') ./ prey
                grazing_zi = replace!(grazing_zi, NaN => 0.0)
                push!(grazing, grazing_zi)
            end   
        end
        return sum(grazing)
    end

    #------- for 1N 4P 3Z 7B 4D
    if nb==7
        np=4
        prey_POM = GrM[2,5]' .* B[:,1]
        gb_POM = g_max .* prey_POM ./ (prey_POM .+ K_g)
        grz_POM = (gb_POM .* Z[:,2] .* GrM[2,5]' ./ prey_POM) 
        push!(grazing, grz_POM)

        prey = GrM[3,6:end]' .* B[:,2:end]
        for i in range(1, nb-1)
            gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
            grz_i = (gb_i .* Z[:,3] .* GrM[3,i+(np+1)]' ./ prey[:,i]) 

            if any(isnan.(grz_i)) 
                return true
            end

            push!(grazing, grz_i)
        end
    
        return grazing
    
    else

        #----- for 1N 2P 2Z 2B 2D
        if nz == 2
            K_g = ds["K_g"][2]
            prey = GrM[2,3:end]' .*B[:,1:end]

            for i in range(1, nb)
                gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
                grz_i = (gb_i .* Z[:,2] .* GrM[2,i+2]' ./ prey[:,i])
                push!(grazing, grz_i)
            end

        #----- for 1N 1P 3Z 2B 2D
        elseif nz == 3
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

    if nb == 13
        loss = zeros(Float64, 89, nb)
        for i in range(1, nb)
            loss[:, i] = mortality[i] .+ grazing[:, i]
        end
    else
        loss = Any[]
        for i in range(1, nb)
            push!(loss, mortality[i] .+ grazing[i])
        end
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

function Rstar(loss, ds, season, nb)

    vmax_ij = ds["vmax_ij"][:]
    Km_ij = ds["Km_ij"][:]
    yield = ds["y_ij"][:]
    temp_mod = get_temp_mod(season)
    II, JJ = get_nonzero_axes(ds["CM"][:])
    RS = Any[]

    if nb == 13
        for j = axes(II, 1)
            push!(RS, Km_ij[II[j],JJ[j]] .* loss[:, j] ./ (yield[II[j],JJ[j]] .* vmax_ij[II[j],JJ[j]] .* temp_mod .- loss[:, j]))
        end
    else
        for j = axes(II, 1)
            push!(RS, Km_ij[II[j],JJ[j]] .* loss[j] ./ (yield[II[j],JJ[j]] .* vmax_ij[II[j],JJ[j]] .* temp_mod .- loss[j]))
        end
    end

    for i in range(1, length(RS))
        RS[i] = check_for_negatives(RS[i])
    end

    return RS

end


#--------------------------------------------------------------------------------------
# R* Phyto on Nutrients
#--------------------------------------------------------------------------------------
function get_rstar_P(P, Z, ds, np, nz, season)
    
    mort_p = p_mortality(P, ds, np)
    grz_p = p_grazing(P, Z, ds, np, nz)
    loss_p = p_loss(mort_p, grz_p, np)
    RstarP_i = RstarP(loss_p, ds, np, season)

    return RstarP_i

end

function p_mortality(P, ds, np)

    m_lp = 1e-1
    m_qp = 1e-1
    mort_p = Any[]
    for i in range(1, np)
        push!(mort_p, (m_lp .+ m_qp .* P[:,i]))
    end

    return mort_p

end

function p_grazing(P, Z, ds, np, nz)

    GrM = ds["GrM"][:]
    grazing = Any[]
    g_max = 1.0
    K_g = ds["K_g"][1]

    if np == 8
        grazing_zi = zeros(Float64, length(P[:,1]), np) 
        for z in range(1, nz)
            if sum(GrM[z,1:np]) > 0 
                prey = GrM[z,1:np]' .* P[:,1:end]
                g = g_max .* prey ./ (prey .+ K_g)
                grazing_zi += (g .* Z[:,z] .* GrM[z,1:np]') ./ prey
                grazing_zi = replace!(grazing_zi, NaN => 0.0)
                push!(grazing, grazing_zi)
            end   
        end

        return sum(grazing)
    
    else

        prey = GrM[1,1:np]' .*P[:,:]
        for i in range(1, np)
            gp_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
            grz_i = (gp_i .* Z[:,1] .* GrM[1,i]') ./ prey[:,i]
            push!(grazing, grz_i)
        end

        return grazing

    end 

end

function p_loss(mortality, grazing, np)

    if np == 8
        loss = zeros(Float64, 89, np)
        for i in range(1, np)
            loss[:, i] = mortality[i] .+ grazing[:, i]
        end
    else
        loss = Any[]
        for i in range(1, np)
            push!(loss, mortality[i] .+ grazing[i])
        end
    end

    return loss

end

function RstarP(loss, ds, np, season)

    umax_ij = ds["umax_ij"][:]
    Kp_ij = ds["Kp_ij"][:]
    temp_mod = get_temp_mod(season)
    RS = Any[]

    if np == 8
        for i in range(1, np)
            push!(RS, Kp_ij[i] .* loss[:, i] ./ (umax_ij[i] .* temp_mod .- loss[:, i]))
        end
    else
        for i in range(1, np)
            push!(RS, (Kp_ij[i] .* loss[i]) ./ (umax_ij[i] .* temp_mod .- loss[i]))
        end
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
function plot_rstar_13B5D(rstar_w, rstar_s, Dw, Ds, bw, bs, ds, N)
    # B1 eats D1 (POM)  |  B2, B6, B10 eat D2  |  B3, B7, B11 eat D3  |  B4, B7, B12 eat D4 | B5, B8, B13 eat D5
    # H = ds["H"][:]
    H = 500
    dz = ds["dz"][:]
    zc = [dz/2:dz:(H-dz/2)]

    # yl=(-400.0, 0)
    lg=:bottomright
    tl=:right
    tfs=22
    # xl=(-0.01, 0.08)
    ab=1.0
    ab_ext=0.3
    ls=5

    p1 = plot(rstar_w[1][1:50], -zc, lw=ls, lc="cyan3", label=L" B1", ylabel="Depth (m)", xrotation=45, xguidefontsize=12, 
    xlabel="", border=:box, legend=lg, xscale=:log10, title="D1", title_loc=tl, titlefontsize=tfs, xlim=(0.008, 0.05))
    plot!(Dw[1:50, 1], -zc, lw=ls, lc="blue3", linestyle=:dot,label=L" D1", alpha=ab)

    p2 = plot(bw[1:50, 1], -zc, lw=ls, lc="cyan3", label=L" B1", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xscale=:log10, xlim=(0.002, 0.02))

    # p3 = plot(rstar_w[2][1:50], -zc, lw=ls, lc="darkorange", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    # border=:box, legend=lg, xscale=:log10, title="D2", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.001, 10))
    p3 = plot(rstar_w[6][1:50], -zc, lw=ls, lc="magenta2", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D2", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.003, 1.12))
    plot!(rstar_w[10][1:50], -zc, lw=ls, lc="purple",  label="")
    plot!(Dw[1:50, 2], -zc, lw=ls, lc="black", linestyle=:dot, alpha=ab, label=L" D2")

    p4 = plot(bw[1:50, 2], -zc, lw=ls, lc="darkorange", label=L" B2", linestyle=:dash, xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab_ext, xlim=(-0.005, 0.034))
    plot!(bw[1:50, 6], -zc, lw=ls, lc="magenta2",  label=L" B6", alpha=ab)
    plot!(bw[1:50, 10], -zc, lw=ls, lc="purple",  label=L" B10", alpha=ab)

    p5 = plot(rstar_w[3][1:50], -zc, lw=ls, lc="gold4", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D3", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.002, 0.014))
    # plot!(rstar_w[7][1:50], -zc, lw=ls, lc="lightcyan4", label="")
    plot!(rstar_w[11][1:50], -zc, lw=ls, lc="darkgrey", label="")
    plot!(Dw[1:50, 3], -zc, lw=ls, lc="maroon", linestyle=:dot, alpha=ab, label=L" D3")

    p6 = plot(bw[1:50, 3], -zc, lw=ls, lc="gold4", label=L" B3", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xlim=(-0.005, 0.032))
    plot!(bw[1:50, 7], -zc, lw=ls, lc="orange4", label=L" B7", linestyle=:dash, alpha=ab_ext)
    plot!(bw[1:50, 11], -zc, lw=ls, lc="darkgrey", label=L" B11", alpha=ab)

    p7 = plot(rstar_w[8][1:50], -zc, lw=ls, lc="seagreen4", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D4", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.001, 0.01))
    # p7 = plot(rstar_w[4][1:50], -zc, lw=ls, lc="coral4", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    # border=:box, legend=lg, xscale=:log10, title="D4", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.001, 0.02))
    # plot!(rstar_w[8][1:50], -zc, lw=ls, lc="seagreen4", label="")
    # plot!(rstar_w[12][1:50], -zc, lw=ls, lc="azure4", label="")
    plot!(Dw[1:50, 4], -zc, lw=ls, lc="navy", linestyle=:dot, alpha=0.6, label=L" D4")

    p8 = plot(bw[1:50, 4], -zc, lw=ls, lc="coral4", linestyle=:dash, label=L" B4", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab_ext, xlim=(-0.005, 0.06))
    plot!(bw[1:50, 8], -zc, lw=ls, lc="seagreen4", label=L" B8", alpha=ab)
    plot!(bw[1:50, 12], -zc, lw=ls, lc="azure4", linestyle=:dash, label=L" B12", alpha=ab_ext)

    p9 = plot(rstar_w[5][1:50], -zc, lw=ls, lc="orchid2", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D5", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.002, 0.01))
    plot!(rstar_w[9][1:50], -zc, lw=ls, lc="darkkhaki", label="")
    # plot!(rstar_w[13][1:50], -zc, lw=ls, lc="turquoise1", label="")
    plot!(Dw[1:50, 5], -zc, lw=ls, lc="brown4", linestyle=:dot, alpha=0.6, label=L" D5")

    p10 = plot(bw[1:50, 5], -zc, lw=ls, lc="orchid2", label=L" B5", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xlim=(-0.005, 0.04))
    plot!(bw[1:50, 9], -zc, lw=ls, lc="darkkhaki", label=L" B9", alpha=ab)
    plot!(bw[1:50, 13], -zc, lw=ls, lc="turquoise1", label=L" B13", linestyle=:dash, alpha=ab_ext)

    #--------------------------------------------------------------------------------------------

    p11 = plot(rstar_s[1][1:50], -zc, lw=ls, lc="cyan3", label=L" B1", ylabel="Depth (m)", xrotation=45, xguidefontsize=12, 
    xlabel=L"R*", border=:box, legend=lg, xscale=:log10, xlim=(0.008, 0.05))
    plot!(Ds[1:50, 1], -zc, lw=ls, lc="blue3", linestyle=:dot,label=L" D1", alpha=ab)

    p12 = plot(bs[1:50, 1], -zc, lw=ls, lc="cyan3", label=L" B1", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xscale=:log10, xlim=(0.002, 0.1))

    # p13 = plot(rstar_s[2][1:50], -zc, lw=ls, lc="darkorange", label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    # border=:box, legend=lg, xscale=:log10, title="D2", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.001, 10))
    p13 = plot(rstar_s[6][1:50], -zc, lw=ls, lc="magenta2", label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), xlim=(0.003, 1.12))
    # plot!(rstar_s[10][1:50], -zc, lw=ls, lc="purple",  label="")
    plot!(Ds[1:50, 2], -zc, lw=ls, lc="black", linestyle=:dot, alpha=0.6, label=L" D2")

    p14 = plot(bs[1:50, 2], -zc, lw=ls, lc="darkorange", label=L" B2", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext, xlim=(-0.01, 0.13))
    plot!(bs[1:50, 6], -zc, lw=ls, lc="magenta2",  label=L" B6", alpha=ab)
    plot!(bs[1:50, 10], -zc, lw=ls, lc="purple",  label=L" B10", linestyle=:dash, alpha=ab_ext)

    # p15 = plot(rstar_s[3][1:50], -zc, lw=ls, lc="gold4", label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    # border=:box, legend=lg, xscale=:log10, title="D3", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.001, 0.1))
    p15 = plot(rstar_s[7][1:50], -zc, lw=ls, lc="lightcyan4", label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), xlim=(0.001, 0.03))
    plot!(rstar_s[11][1:50], -zc, lw=ls, lc="orange4", label="")
    plot!(Ds[1:50, 3], -zc, lw=ls, lc="maroon", linestyle=:dot, alpha=0.6, label=L" D3")

    p16 = plot(bs[1:50, 3], -zc, lw=ls, lc="gold4", label=L" B3", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext, xlim=(-0.005, 0.09))
    plot!(bs[1:50, 7], -zc, lw=ls, lc="orange4", label=L" B7", alpha=ab)
    plot!(bs[1:50, 11], -zc, lw=ls, lc="darkgrey", label=L" B11", alpha=ab)

    p17 = plot(rstar_s[4][1:50], -zc, lw=ls, lc="coral4", label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), xlim=(0.001, 0.02))
    plot!(rstar_s[8][1:50], -zc, lw=ls, lc="seagreen4", label="")
    plot!(rstar_s[12][1:50], -zc, lw=ls, lc="azure4", label="")
    plot!(Ds[1:50, 4], -zc, lw=ls, lc="navy", linestyle=:dot, alpha=0.6, label=L" D4")

    p18 = plot(bs[1:50, 4], -zc, lw=ls, lc="coral4", label=L" B4", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xlim=(-0.005, 0.08))
    plot!(bs[1:50, 8], -zc, lw=ls, lc="seagreen4", label=L" B8", alpha=ab)
    plot!(bs[1:50, 12], -zc, lw=ls, lc="azure4", label=L" B12", alpha=ab)

    p19 = plot(rstar_s[5][1:50], -zc, lw=ls, lc="orchid2", label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), xlim=(0.002, 0.02))
    plot!(rstar_s[9][1:50], -zc, lw=ls, lc="darkkhaki", label="")
    # plot!(rstar_s[13][1:50], -zc, lw=ls, lc="turquoise1", label="")
    plot!(Ds[1:50, 5], -zc, lw=ls, lc="brown4", linestyle=:dot, alpha=ab, label=L" D5")

    p20 = plot(bs[1:50, 5], -zc, lw=ls, lc="orchid2", label=L" B5", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xlim=(-0.01, 0.11))
    plot!(bs[1:50, 9], -zc, lw=ls, lc="darkkhaki", label=L" B9", alpha=ab)
    plot!(bs[1:50, 13], -zc, lw=ls, lc="turquoise1", linestyle=:dash, label=L" B13", alpha=ab_ext)
    
    winter = plot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
    fg_legend = :transparent,
    layout = (1,10),
    size=(1500,450),
    # xlabel = "R*",
    # plot_title="Winter", 
    # plot_titlefontsize = 20,
    # titlefontsize=tfs, titlelocation=:center, 
    )

    summer = plot(p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,
    fg_legend = :transparent,
    layout = (1,10),
    size=(1500,450),
    # xlabel = "R*",
    # plot_title="Summer", 
    # plot_titlefontsize = 20,
    # titlefontsize=tfs, titlelocation=:center, 
    )
    
    combined = plot(winter, summer, 
        fg_legend = :transparent,
        layout = (2,1),
        size=(1500,900),
        plot_xaxis="Depth (m)"
        # xlabel = "R*",
        # plot_title="Winter", 
        # titlefontsize=tfs, titlelocation=:center, 
    )

    savefig(combined,"/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rstarB_13B5D.png")
    
    return combined

end

function plot_rstar_7B4D(rstar_w, rstar_s, Dw, Ds, Bw, Bs, ds)
    # B1 eats D1 (POM)  |  B2 & B7 eat D2  |  B3 & B6 eat D3  |  B4 & B5 eat D4

    H = 500
    dz = 10
    zc = [dz/2:dz:(H-dz/2)]
    lg=:bottomright
    tfs=22
    tl=:right
    ab=1.0
    ab_ext=0.3
    ls=5

    #--------------------------------------------------------------------------------------------
    # Winter
    #--------------------------------------------------------------------------------------------
    p1 = plot(rstar_w[1][1:50], -zc, lw=ls, lc="red3", label=L" B1", ylabel="Depth (m)", xrotation=45, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D1", title_loc=tl, titlefontsize=tfs)
    plot!(Dw[1:50, 1], -zc, lw=ls, lc="blue3", linestyle=:dot, label=L" D1", alpha=ab)

    p2 = plot(Bw[1:50, 1], -zc, lw=ls, lc="cyan3", label=L" B1", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xscale=:log10)

    p3 = plot(rstar_w[6][1:50], -zc, lw=ls, lc="magenta2", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), title="D3", title_loc=tl, titlefontsize=tfs)
    # plot!(rstar_w[3][1:50], -zc, lw=ls, lc="darkorange",  label="")
    plot!(Dw[1:50, 3], -zc, lw=ls, lc="black", linestyle=:dot, alpha=ab, label=L" D3")

    p4 = plot(Bw[1:50, 3], -zc, lw=ls, lc="darkorange", label=L" B3", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext)
    plot!(Bw[1:50, 6], -zc, lw=ls, lc="magenta2",  label=L" B6", alpha=ab)

    p5 = plot(rstar_w[5][1:50], -zc, lw=ls, lc="lightcyan4", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), title="D4", title_loc=tl, titlefontsize=tfs)
    # plot!(rstar_w[4][1:50], -zc, lw=ls, lc="gold4", label="")
    plot!(Dw[1:50, 4], -zc, lw=ls, lc="purple", linestyle=:dot, alpha=ab, label=L" D4")

    p6 = plot(Bw[1:50, 4], -zc, lw=ls, lc="gold4", label=L" B4", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext)
    plot!(Bw[1:50, 5], -zc, lw=ls, lc="lightcyan4", label=L" B5", alpha=ab)

    #--------------------------------------------------------------------------------------------
    # Summer
    #--------------------------------------------------------------------------------------------
    p7 = plot(rstar_s[1][1:50], -zc, lw=ls, lc="red3", label=L" B1", ylabel="Depth (m)", 
    xrotation=45, xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10)
    plot!(Ds[1:50, 1], -zc, lw=ls, lc="blue3", linestyle=:dot,label=L" D1", alpha=ab)

    p8 = plot(Bs[1:50, 1], -zc, lw=ls, lc="cyan3", label=L" B1", xrotation=45, xguidefontsize=12, 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xscale=:log10, xlabel=L"Biomass")

    p9 = plot(rstar_s[6][1:50], -zc, lw=ls, lc="magenta2", label="", xrotation=45, xguidefontsize=12, 
    xlabel=L"R*", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), alpha=ab)
    # plot!(rstar_s[3][1:50], -zc, lw=ls, lc="darkorange", label="")
    plot!(Ds[1:50, 3], -zc, lw=ls, lc="black", linestyle=:dot, alpha=ab, label=L" D3")

    p10 = plot(Bs[1:50, 3], -zc, lw=ls, lc="darkorange", label=L" B3", xrotation=45, xguidefontsize=12, 
    xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext)
    plot!(Bs[1:50, 6], -zc, lw=ls, lc="magenta2", label=L" B6", alpha=ab)

    p11 = plot(rstar_s[4][1:50], -zc, lw=ls, lc="gold4", label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    # plot!(rstar_s[5][1:50], -zc, lw=ls, lc="lightcyan4", label="")
    plot!(Ds[1:50, 4], -zc, lw=ls, lc="purple", linestyle=:dot, alpha=0.6, label=L" D4")

    p12 = plot(Bs[1:50, 4], -zc, lw=ls, lc="gold4", label=L" B4", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""))
    plot!(Bs[1:50, 5], -zc, lw=ls, lc="lightcyan4", label=L" B5", linestyle=:dash, alpha=ab_ext)
    
    
    combined = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
        fg_legend = :transparent,
        layout = (2,6),
        size=(900,700),
    )

    savefig(combined, "/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rstarB_7B4D.png")

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
    dz = 10
    zc = [dz/2:dz:(H-dz/2)]
    lg=:bottomright
    tfs=12
    al=0.7
    al_ext=0.3
    linw=5

    if np == 1
        p1 = plot(rstar_w[1][1:20], -zc, lw=linw, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, title="", alpha=al)
        plot!(Nw[1:20, 1], -zc, lw=linw, lc="blue3", linestyle=:dot, label=L" N", legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=linw, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))

        p3 = plot(rstar_s[1][1:20], -zc, lw=linw, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, title="", alpha=al)
        plot!(Ns[1:20, 1], -zc, lw=linw, lc="blue3", linestyle=:dot,label=L" N", legend=:bottom)

        p4 = plot(Ps[1:20, 1], -zc, lw=linw, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))
    end
    
    if np == 2
        p1 = plot(rstar_w[1][1:20], -zc, lw=linw, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(rstar_w[2][1:20], -zc, lw=linw, lc="darkgreen",  label="")
        plot!(Nw[1:20, 1], -zc, lw=linw, lc="blue3", linestyle=:dot, label=L" N", alpha=0.4, legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=linw, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))
        plot!(Pw[1:20, 2], -zc, lw=linw, lc="darkgreen",  label=L" P2")

        p3 = plot(rstar_s[1][1:20], -zc, lw=linw, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(rstar_s[2][1:20], -zc, lw=linw, lc="darkgreen",  label="")
        plot!(Ns[1:20, 1], -zc, lw=linw, lc="blue3", linestyle=:dot,label=L" N", alpha=0.4, legend=:bottom)

        p4 = plot(Ps[1:20, 1], -zc, lw=linw, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))
        plot!(Ps[1:20, 2], -zc, lw=linw, lc="darkgreen",  label=L" P2")
    end 

    if np == 4
        p1 = plot(rstar_w[2][1:20], -zc, lw=linw, lc="darkgreen", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10)
        # plot!(rstar_w[1][1:20], -zc, lw=linw, lc="olivedrab",  label="")
        # plot!(rstar_w[3][1:20], -zc, lw=linw, lc="red4",  label="")
        plot!(rstar_w[4][1:20], -zc, lw=linw, lc="cyan4",  label="")
        plot!(Nw[1:20, 1], -zc, lw=linw, lc="blue3", linestyle=:dot, label=L" N", alpha=al, legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=linw, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=al_ext)
        plot!(Pw[1:20, 2], -zc, lw=linw, lc="darkgreen",  label=L" P2")
        plot!(Pw[1:20, 3], -zc, lw=linw, lc="red4",  label=L" P3", linestyle=:dash, alpha=al_ext)
        plot!(Pw[1:20, 4], -zc, lw=linw, lc="cyan4",  label=L" P4")

        p3 = plot(rstar_s[1][1:20], -zc, lw=linw, lc="olivedrab3", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        # plot!(rstar_s[2][1:20], -zc, lw=linw, lc="darkgreen",  label="")
        plot!(rstar_s[3][1:20], -zc, lw=linw, lc="red4",  label="")
        # plot!(rstar_s[4][1:20], -zc, lw=linw, lc="cyan4",  label="")
        plot!(Ns[1:20, 1], -zc, lw=linw, lc="blue3", linestyle=:dot,label=L" N", alpha=0.4, legend=:bottom)

        p4 = plot(Ps[1:20, 1], -zc, lw=linw, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), xlims=(-0.01, 0.6))
        plot!(Ps[1:20, 2], -zc, lw=linw, lc="darkgreen",  label=L" P2", linestyle=:dash, alpha=al_ext)
        plot!(Ps[1:20, 3], -zc, lw=linw, lc="red4",  label=L" P3")
        plot!(Ps[1:20, 4], -zc, lw=linw, lc="cyan4",  label=L" P4", linestyle=:dash, alpha=al_ext) 
    end

    if np == 8
        p1 = plot(rstar_w[2][1:20], -zc, lw=linw, lc="darkgreen", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, alpha=al)
        # plot!(rstar_w[1][1:20], -zc, lw=linw, lc="dolivedrab",  label="", alpha=al)
        # plot!(rstar_w[3][1:20], -zc, lw=linw, lc="red4",  label="", alpha=al)
        # plot!(rstar_w[4][1:20], -zc, lw=linw, lc="cyan4",  label="", alpha=al)
        # plot!(rstar_w[5][1:20], -zc, lw=linw, lc="purple",  label="", alpha=al)
        plot!(rstar_w[6][1:20], -zc, lw=linw, lc="black",  label="", alpha=al)
        # plot!(rstar_w[7][1:20], -zc, lw=linw, lc="hotpink2",  label="", alpha=al)
        # plot!(rstar_w[8][1:20], -zc, lw=linw, lc="wheat2",  label="", alpha=al)
        plot!(Nw[1:20, 1], -zc, lw=linw, lc="blue3", linestyle=:dot, label=L" N", legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=linw, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, linestyle=:dash, alpha=al_ext,  yformatter=Returns(""), xlims=(-0.01, 0.15))
        plot!(Pw[1:20, 2], -zc, lw=linw, lc="darkgreen", label=L" P2")
        plot!(Pw[1:20, 3], -zc, lw=linw, lc="red4", linestyle=:dash, alpha=al_ext,  label=L" P3")
        plot!(Pw[1:20, 4], -zc, lw=linw, lc="cyan4", linestyle=:dash, alpha=al_ext,   label=L" P4")
        plot!(Pw[1:20, 5], -zc, lw=linw, lc="purple", linestyle=:dash, alpha=al_ext,   label=L" P5")
        plot!(Pw[1:20, 6], -zc, lw=linw, lc="black",  label=L" P6")
        plot!(Pw[1:20, 7], -zc, lw=linw, lc="hotpink2", linestyle=:dash, alpha=al_ext,   label=L" P7")
        plot!(Pw[1:20, 8], -zc, lw=linw, lc="wheat2", linestyle=:dash, alpha=al_ext,   label=L" P8")

        p3 = plot(rstar_s[4][1:20], -zc, lw=linw, lc="cyan4", label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, alpha=al)
        # plot!(rstar_s[1][1:20], -zc, lw=linw, lc="olivedrab3",  label="", alpha=al)
        # plot!(rstar_s[2][1:20], -zc, lw=linw, lc="darkgreen",  label="", alpha=al)
        # plot!(rstar_s[3][1:20], -zc, lw=linw, lc="red4",  label="", alpha=al)
        # plot!(rstar_s[4][1:20], -zc, lw=linw, lc="cyan4",  label="", alpha=al)
        plot!(rstar_s[5][1:20], -zc, lw=linw, lc="purple",  label="", alpha=al)
        # plot!(rstar_s[6][1:20], -zc, lw=linw, lc="black",  label="", alpha=al)
        # plot!(rstar_s[7][1:20], -zc, lw=linw, lc="hotpink2",  label="", alpha=al)
        # plot!(rstar_s[8][1:20], -zc, lw=linw, lc="wheat2",  label="", alpha=al)
        plot!(Ns[1:20, 1], -zc, lw=linw, lc="blue3", linestyle=:dot,label=L" N", legend=:bottom)

        p4 = plot(Ps[1:20, 1], -zc, lw=linw, lc="olivedrab3", label=L" P1", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=al_ext)
        plot!(Ps[1:20, 2], -zc, lw=linw, lc="darkgreen",  label=L" P2", linestyle=:dash, alpha=al_ext)
        plot!(Ps[1:20, 3], -zc, lw=linw, lc="red4",  label=L" P3", linestyle=:dash, alpha=al_ext)
        plot!(Ps[1:20, 4], -zc, lw=linw, lc="cyan4",  label=L" P4") 
        plot!(Ps[1:20, 5], -zc, lw=linw, lc="purple",  label=L" P5")
        plot!(Ps[1:20, 6], -zc, lw=linw, lc="black",  label=L" P6", linestyle=:dash, alpha=al_ext)
        plot!(Ps[1:20, 7], -zc, lw=linw, lc="hotpink2",  label=L" P7", linestyle=:dash, alpha=al_ext) 
        plot!(Ps[1:20, 8], -zc, lw=linw, lc="wheat2",  label=L" P8", linestyle=:dash, alpha=al_ext) 

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
        RS[i] = ifelse(RS[i] > 5, NaN, RS[i])
    end

    return RS

end

function get_zc(ds, H=890)

    dz = 10
    zc = [dz/2:dz:(H-dz/2)]

    return zc
end



#################################  RUN  #####################################

# -----------------------------------------------------------------------
# 1N 8P 6Z 13B 5D
# -----------------------------------------------------------------------
# winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_230915_22:19_8P6Z13B5D.nc")
# summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Su100y_230916_11:25_8P6Z13B5D.nc")
# CMd = display(winter["CM"][:])
# GrMd = display(winter["GrM"][:])

# # Get endpoints and set all biomass values below cut-off (10^-6) equal to cut-off
# Nw, Pw, Zw, Bw, Dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
# Ns, Ps, Zs, Bs, Ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

# Bww = cut_off(Bw, 13)
# Bss = cut_off(Bs, 13)
# Pww = cut_off(Pw, 8)
# Pss = cut_off(Ps, 8)

# # Calculate Rstar B
# rstarB_w = get_rstar_B(Bww, Zw, winter, 13, 6, "Win")
# rstarB_s = get_rstar_B(Bss, Zs, summer, 13, 6, "Sum")
# rstarP_w = get_rstar_P(Pww, Zw, winter, 8, 6, "Win")
# rstarP_s = get_rstar_P(Pss, Zs, summer, 8, 6, "Sum")

# Bwe = extinct(Bww, 13)
# Bse = extinct(Bss, 13)
# Pwe = extinct(Pww, 8)
# Pse = extinct(Pss, 8)

# plot_rstar_13B5D(rstarB_w, rstarB_s, Dw, Ds, Bwe, Bse, winter, 13)
# plot_rstar_P(rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, winter, 8)

# # -----------------------------------------------------------------------
# # 1N 4P 3Z 7B 4D
# # -----------------------------------------------------------------------
winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230827_13:45_4P3Z7B4D_ep.nc")
summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Su100y_230827_17:10_4P3Z7B4D_ep.nc")
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
rstarP_w = get_rstar_P(Pww, Zw, winter, 4, 3, "Win")
rstarP_s = get_rstar_P(Pss, Zs, summer, 4, 3, "Sum")

Bwe = extinct(Bww, 7)
Bse = extinct(Bss, 7)
Pwe = extinct(Pww, 4)
Pse = extinct(Pss, 4)

plot_rstar_7B4D(rstarB_ijw, rstarB_ijs, Dw, Ds, Bwe, Bse, winter)
plot_rstar_P(rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, winter,  4)
println("done")


# # -----------------------------------------------------------------------
# # 1N 2P 2Z 2B 2D
# # -----------------------------------------------------------------------
# winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230905_2005.nc")
# summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230905_2301.nc")
# CMd = display(winter["CM"][:])
# GrMd = display(winter["GrM"][:])

# # Get endpoints 
# Nw, Pw, Zw, Bw, Dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
# Ns, Ps, Zs, Bs, Ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

# # set all biomass values below cut-off (10^-6) equal to cut-off
# Bww = cut_off(Bw, 2)
# Bss = cut_off(Bs, 2)
# Pww = cut_off(Pw, 2)
# Pss = cut_off(Ps, 2)

# # Calculate Rstar 
# rstarB_ijw = get_rstar_B(Bww, Zw, winter, 2, 2, "Win")
# rstarB_ijs = get_rstar_B(Bss, Zs, summer, 2, 2, "Sum")
# rstarP_w = get_rstar_P(Pww, Zw, winter, 2, "Win")
# rstarP_s = get_rstar_P(Pss, Zs, summer, 2, "Sum")

# # set low biomass to zero
# Bwe = extinct(Bww, 2)
# Bse = extinct(Bss, 2)
# Pwe = extinct(Pww, 2)
# Pse = extinct(Pss, 2)


# # Plot
# plot_rstar_2B2D(rstarB_ijw, rstarB_ijs, Dw, Ds, Bwe, Bse, winter, 2)
# plot_rstar_P(rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, winter,  2)
# println("done")


# # -----------------------------------------------------------------------
# # 1N 1P 3Z 2B 2D
# # -----------------------------------------------------------------------
# winter = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230906_1351.nc")
# summer = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/out_100y_20230906_1413.nc")
# CMd = display(winter["CM"][:])
# GrMd = display(winter["GrM"][:])

# # Get endpoints
# Nw, Pw, Zw, Bw, Dw = get_endpoints(winter, ["n", "p", "z", "b", "d"])
# Ns, Ps, Zs, Bs, Ds = get_endpoints(summer, ["n", "p", "z", "b", "d"])

# # set all biomass values below cut-off (10^-6) equal to cut-off
# Bww = cut_off(Bw, 2)
# Bss = cut_off(Bs, 2)
# Pww = cut_off(Pw, 1)
# Pss = cut_off(Ps, 1)

# # Calculate Rstar 
# rstarB_ijw = get_rstar_B(Bww, Zw, winter, 2, 3, "Win")
# rstarB_ijs = get_rstar_B(Bss, Zs, summer, 2, 3, "Sum")
# rstarPw = get_rstar_P(Pww, Zw, winter, 1, "Win")
# rstarPs = get_rstar_P(Pss, Zs, summer, 1, "Sum")

# Bwe = extinct(Bww, 2)
# Bse = extinct(Bss, 2)
# Pwe = extinct(Pww, 1)
# Pse = extinct(Pss, 1)

# # Plot
# plot_rstar_2B2D(rstarB_ijw, rstarB_ijs, Dw, Ds, Bwe, Bse, winter, 1)
# plot_rstar_P(rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, winter,  1)
# println("done")



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