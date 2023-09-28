using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames, CSV
using SparseArrays, LinearAlgebra


#--------------------------------------------------------------------------------------
# R* Bacteria on Detritus
#--------------------------------------------------------------------------------------
# function get_rstar_B(B, Z, ds, nb, nz, season)
    
#     mort_b = b_mortality(B, ds, nb)
#     grz_b = b_grazing(B, Z, ds, nb, nz)
#     loss_b = b_loss(mort_b, grz_b, nb)
#     RstarB_ij = Rstar(loss_b, ds, season, nb)

#     return RstarB_ij

# end


# function b_mortality(B, ds, nb)

#     mort_b = Any[]
#     for i in range(1, nb)
#         push!(mort_b, (ds["m_lb"][i] .+ ds["m_qb"][i] .* B[:,i]))
#     end

#     return mort_b

# end

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

# function b_grazing(B, Z, ds, nb, nz)
#     #NOTE this all needs to be generalized / scaled !! 
#     GrM = ds["GrM"][:]
#     grazing = Any[]
#     g_max = 1.0
#     K_g = 1.0

#     #------- for 1N 8P 6Z 13B 5D
#     if nb == 13
#         # NOTE length(B[:,1]) = prms.ngrid
#         grazing_zi = zeros(Float64, length(B[:,1]), nb) 
#         np = 8
#         for z in range(1, nz)
#             if sum(GrM[z,np+1:end]) > 0 
#                 prey = GrM[z,np+1:end]' .*B[:,1:end]
#                 g = g_max .* prey ./ (prey .+ K_g)
#                 grazing_zi += (g .* Z[:,z] .* GrM[z,np+1:end]') ./ prey
#                 grazing_zi = replace!(grazing_zi, NaN => 0.0)
#                 push!(grazing, grazing_zi)
#             end   
#         end
#         return sum(grazing)
#     end

#     #------- for 1N 4P 3Z 7B 4D
#     if nb==7
#         np=4
#         prey_POM = GrM[2,5]' .* B[:,1]
#         gb_POM = g_max .* prey_POM ./ (prey_POM .+ K_g)
#         grz_POM = (gb_POM .* Z[:,2] .* GrM[2,5]' ./ prey_POM) 
#         push!(grazing, grz_POM)

#         prey = GrM[3,6:end]' .* B[:,2:end]
#         for i in range(1, nb-1)
#             gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
#             grz_i = (gb_i .* Z[:,3] .* GrM[3,i+(np+1)]' ./ prey[:,i]) 

#             if any(isnan.(grz_i)) 
#                 return true
#             end

#             push!(grazing, grz_i)
#         end
    
#         return grazing
    
#     else

#         #----- for 1N 2P 2Z 2B 2D
#         if nz == 2
#             K_g = ds["K_g"][2]
#             prey = GrM[2,3:end]' .*B[:,1:end]

#             for i in range(1, nb)
#                 gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
#                 grz_i = (gb_i .* Z[:,2] .* GrM[2,i+2]' ./ prey[:,i])
#                 push!(grazing, grz_i)
#             end

#         #----- for 1N 1P 3Z 2B 2D
#         elseif nz == 3
#             K_g = [ds["K_g"][2], ds["K_g"][3]]
#             prey = 1.0 .*B[:,1:2]

#             for i in range(1, nb)
#                 gb_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g[i])
#                 grz_i = (gb_i .* Z[:,i+1] .* 1.0 ./ prey[:,i]) 
#                 push!(grazing, grz_i)
#             end
#         else
#         end
#     end

#     return grazing

# end

# function b_loss(mortality, grazing, nb)

#     if nb == 13
#         loss = zeros(Float64, 89, nb)
#         for i in range(1, nb)
#             loss[:, i] = mortality[i] .+ grazing[:, i]
#         end
#     else
#         loss = Any[]
#         for i in range(1, nb)
#             push!(loss, mortality[i] .+ grazing[i])
#         end
#     end
#     return loss
# end

# function get_temp_mod(season)
#     #fit to SPOT data (approx 20 to 4, approx 16 to 4)
#     if season == "Win"
#         temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/win_temp_mod.csv", DataFrame)
#     else
#         temp_mod = CSV.read("/home/lee/Dropbox/Development/NPZBD_1D/data/temp_mod/sum_temp_mod.csv", DataFrame)
#     end

#     return Matrix(temp_mod)
# end

# function Rstar(loss, ds, season, nb)

#     vmax_ij = ds["vmax_ij"][:]
#     Km_ij = ds["Km_ij"][:]
#     yield = ds["y_ij"][:]
#     temp_mod = get_temp_mod(season)
#     II, JJ = get_nonzero_axes(ds["CM"][:])
#     RS = Any[]

#     if nb == 13
#         for j = axes(II, 1)
#             push!(RS, Km_ij[II[j],JJ[j]] .* loss[:, j] ./ (yield[II[j],JJ[j]] .* vmax_ij[II[j],JJ[j]] .* temp_mod .- loss[:, j]))
#         end
#     else
#         for j = axes(II, 1)
#             push!(RS, Km_ij[II[j],JJ[j]] .* loss[j] ./ (yield[II[j],JJ[j]] .* vmax_ij[II[j],JJ[j]] .* temp_mod .- loss[j]))
#         end
#     end

#     for i in range(1, length(RS))
#         RS[i] = check_for_negatives(RS[i])
#     end

#     return RS

# end


#--------------------------------------------------------------------------------------
# R* Phyto on Nutrients
#--------------------------------------------------------------------------------------
# function get_rstar_P(P, Z, ds, np, nz, season)
    
#     mort_p = p_mortality(P, ds, np)
#     grz_p = p_grazing(P, Z, ds, np, nz)
#     loss_p = p_loss(mort_p, grz_p, np)
#     RstarP_i = RstarP(loss_p, ds, np, season)

#     return RstarP_i

# end

# function p_mortality(P, ds, np)

#     m_lp = 1e-1
#     m_qp = 1e-1
#     mort_p = Any[]
#     for i in range(1, np)
#         push!(mort_p, (m_lp .+ m_qp .* P[:,i]))
#     end

#     return mort_p

# end

# function p_grazing(P, Z, ds, np, nz)

#     GrM = ds["GrM"][:]
#     grazing = Any[]
#     g_max = 1.0
#     K_g = ds["K_g"][1]

#     if np == 8
#         grazing_zi = zeros(Float64, length(P[:,1]), np) 
#         for z in range(1, nz)
#             if sum(GrM[z,1:np]) > 0 
#                 prey = GrM[z,1:np]' .* P[:,1:end]
#                 g = g_max .* prey ./ (prey .+ K_g)
#                 grazing_zi += (g .* Z[:,z] .* GrM[z,1:np]') ./ prey
#                 grazing_zi = replace!(grazing_zi, NaN => 0.0)
#                 push!(grazing, grazing_zi)
#             end   
#         end

#         return sum(grazing)
    
#     else

#         prey = GrM[1,1:np]' .*P[:,:]
#         for i in range(1, np)
#             gp_i = g_max .* prey[:,i] ./ (prey[:,i] .+ K_g)
#             grz_i = (gp_i .* Z[:,1] .* GrM[1,i]') ./ prey[:,i]
#             push!(grazing, grz_i)
#         end

#         return grazing

#     end 

# end

# function p_loss(mortality, grazing, np)

#     if np == 8
#         loss = zeros(Float64, 89, np)
#         for i in range(1, np)
#             loss[:, i] = mortality[i] .+ grazing[:, i]
#         end
#     else
#         loss = Any[]
#         for i in range(1, np)
#             push!(loss, mortality[i] .+ grazing[i])
#         end
#     end

#     return loss

# end

# function RstarP(loss, ds, np, season)

#     umax_ij = ds["umax_ij"][:]
#     Kp_ij = ds["Kp_ij"][:]
#     temp_mod = get_temp_mod(season)
#     RS = Any[]

#     if np == 8
#         for i in range(1, np)
#             push!(RS, Kp_ij[i] .* loss[:, i] ./ (umax_ij[i] .* temp_mod .- loss[:, i]))
#         end
#     else
#         for i in range(1, np)
#             push!(RS, (Kp_ij[i] .* loss[i]) ./ (umax_ij[i] .* temp_mod .- loss[i]))
#         end
#     end

#     for i in range(1, length(RS))
#         RS[i] = check_for_negatives(RS[i])
#     end

#     return RS

# end

#--------------------------------------------------------------------------------------
# R* Z on B & P
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
# Plotting
#--------------------------------------------------------------------------------------
bc= ["cyan3", "darkorange", "indigo", "coral4", "lightcyan4", "magenta2", "orange4", "seagreen4",
"darkkhaki", "purple", "crimson",  "azure4", "turquoise1"]
dc= ["blue3", "black", "maroon", "navy", "brown4"]
pc = ["olivedrab3", "darkgreen","red4", "cyan4", "purple", "black", "hotpink2", "wheat2" ]
nc = ["blue2"]
ab=0.8
ab_ext=0.8
ls=5
lfs=9

function plot_rstar_13B5D(fsaven, rstar_w, rstar_s, Dw, Ds, bw, bs)
    # B1 eats D1 (POM)  |  B2, B6, B10 eat D2  |  B3, B7, B11 eat D3  |  B4, B7, B12 eat D4 | B5, B8, B13 eat D5
    # H = ds["H"][:]

    H = 500
    dz = 10
    zc = [dz/2:dz:(H-dz/2)]
    # yl=(-400.0, 0)
    lg=:bottomright
    tl=:right
    tfs=22
    lfs=9

    p1 = plot(rstar_w[1][1:50], -zc, lw=ls, lc=bc[1], label=" B1", ylabel="Depth (m)", xrotation=45, xguidefontsize=12, 
    xlabel="", border=:box, legend=lg, xscale=:log10, title="D1", title_loc=tl, titlefontsize=tfs, legendfontsize=lfs)
    plot!(Dw[1:50, 1], -zc, lw=ls, lc=dc[1], linestyle=:dot,label=" D1", alpha=ab, legendfontsize=lfs)

    p2 = plot(bw[1:50, 1], -zc, lw=ls, lc=bc[1], label=" B1", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xscale=:log10, legendfontsize=lfs)

    # p3 = plot(rstar_w[2][1:50], -zc, lw=ls, lc=bc[2], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    # border=:box, legend=lg, xscale=:log10, title="D2", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    p3 = plot(rstar_w[6][1:50], -zc, lw=ls, lc=bc[6], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D2", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    plot!(rstar_w[10][1:50], -zc, lw=ls, lc=bc[10],  label="")
    plot!(Dw[1:50, 2], -zc, lw=ls, lc=dc[2], linestyle=:dot, alpha=ab, label=" D2", legendfontsize=lfs)

    p4 = plot(bw[1:50, 2], -zc, lw=ls, lc=bc[2], label=" B2", linestyle=:dash, xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab_ext, legendfontsize=lfs)
    plot!(bw[1:50, 6], -zc, lw=ls, lc=bc[6],  label=" B6", alpha=ab, legendfontsize=lfs)
    plot!(bw[1:50, 10], -zc, lw=ls, lc=bc[10],  label=" B10", alpha=ab, legendfontsize=lfs)

    p5 = plot(rstar_w[3][1:50], -zc, lw=ls, lc=bc[3], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D3", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    # plot!(rstar_w[7][1:50], -zc, lw=ls, lc=bc[7], label="")
    plot!(rstar_w[11][1:50], -zc, lw=ls, lc=bc[11], label="")
    plot!(Dw[1:50, 3], -zc, lw=ls, lc=dc[3], linestyle=:dot, alpha=ab, label=" D3", legendfontsize=lfs)

    p6 = plot(bw[1:50, 3], -zc, lw=ls, lc=bc[3], label=" B3", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab,legendfontsize=lfs)
    plot!(bw[1:50, 7], -zc, lw=ls, lc=bc[7], label=" B7", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
    plot!(bw[1:50, 11], -zc, lw=ls, lc=bc[11], label=" B11", alpha=ab, legendfontsize=lfs)

    p7 = plot(rstar_w[8][1:50], -zc, lw=ls, lc=bc[8], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D4", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    # p7 = plot(rstar_w[4][1:50], -zc, lw=ls, lc="coral4", label="", xrotation=45, xguidefontsize=12, xlabel="", 
    # border=:box, legend=lg, xscale=:log10, title="D4", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.001, 0.02))
    # plot!(rstar_w[8][1:50], -zc, lw=ls, lc=bc[8], label="")
    # plot!(rstar_w[12][1:50], -zc, lw=ls, lc="azure4", label="")
    plot!(Dw[1:50, 4], -zc, lw=ls, lc=dc[4], linestyle=:dot, alpha=0.6, label=" D4", legendfontsize=lfs)

    p8 = plot(bw[1:50, 4], -zc, lw=ls, lc=bc[4], linestyle=:dash, label=" B4", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab_ext,  legendfontsize=lfs)
    plot!(bw[1:50, 8], -zc, lw=ls, lc=bc[8], label=" B8", alpha=ab, legendfontsize=lfs)
    plot!(bw[1:50, 12], -zc, lw=ls, lc=bc[12], linestyle=:dash, label=" B12", alpha=ab_ext, legendfontsize=lfs)

    p9 = plot(rstar_w[5][1:50], -zc, lw=ls, lc=bc[5], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D5", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    plot!(rstar_w[9][1:50], -zc, lw=ls, lc=bc[9], label="")
    # plot!(rstar_w[13][1:50], -zc, lw=ls, lc="turquoise1", label="")
    plot!(Dw[1:50, 5], -zc, lw=ls, lc=dc[5], linestyle=:dot, alpha=0.6, label=" D5", legendfontsize=lfs)

    p10 = plot(bw[1:50, 5], -zc, lw=ls, lc=bc[5], label=" B5", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, legendfontsize=lfs)
    plot!(bw[1:50, 9], -zc, lw=ls, lc=bc[9], label=" B9", alpha=ab, legendfontsize=lfs)
    plot!(bw[1:50, 13], -zc, lw=ls, lc=bc[13], label=" B13", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)

    #--------------------------------------------------------------------------------------------

    p11 = plot(rstar_s[1][1:50], -zc, lw=ls, lc=bc[1], label=" B1", ylabel="Depth (m)", xrotation=45, xguidefontsize=12, 
    xlabel=L"R*", border=:box, legend=lg, xscale=:log10)
    plot!(Ds[1:50, 1], -zc, lw=ls, lc=dc[1], linestyle=:dot,label=" D1", alpha=ab, legendfontsize=lfs)

    p12 = plot(bs[1:50, 1], -zc, lw=ls, lc=bc[1], label=" B1", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xscale=:log10, legendfontsize=lfs)

    # p13 = plot(rstar_s[2][1:50], -zc, lw=ls, lc=bc[2], label="", xrotation=45, xguidefontsize=12, xlabel="R*", 
    # border=:box, legend=lg, xscale=:log10, title="D2", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs, xlim=(0.001, 10))
    p13 = plot(rstar_s[6][1:50], -zc, lw=ls, lc=bc[6], label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    # plot!(rstar_s[10][1:50], -zc, lw=ls, lc=bc[10],  label="")
    plot!(Ds[1:50, 2], -zc, lw=ls, lc=dc[2], linestyle=:dot, alpha=0.6, label=" D2", legendfontsize=lfs)

    p14 = plot(bs[1:50, 2], -zc, lw=ls, lc=bc[2], label=" B2", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
    plot!(bs[1:50, 6], -zc, lw=ls, lc=bc[6],  label=" B6", alpha=ab, legendfontsize=lfs)
    plot!(bs[1:50, 10], -zc, lw=ls, lc=bc[10],  label=" B10", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)

    # p15 = plot(rstar_s[3][1:50], -zc, lw=ls, lc=bc[3], label="", xrotation=45, xguidefontsize=12, xlabel="R*", 
    # border=:box, legend=lg, xscale=:log10, title="D3", title_loc=tl, yformatter=Returns(""), titlefontsize=tfs)
    p15 = plot(rstar_s[7][1:50], -zc, lw=ls, lc=bc[7], label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    plot!(rstar_s[11][1:50], -zc, lw=ls, lc=bc[11], label="")
    plot!(Ds[1:50, 3], -zc, lw=ls, lc=dc[3], linestyle=:dot, alpha=0.6, label=" D3", legendfontsize=lfs)

    p16 = plot(bs[1:50, 3], -zc, lw=ls, lc=bc[3], label=" B3", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
    plot!(bs[1:50, 7], -zc, lw=ls, lc=bc[7], label=" B7", alpha=ab, legendfontsize=lfs)
    plot!(bs[1:50, 11], -zc, lw=ls, lc=bc[11], label=" B11", alpha=ab, legendfontsize=lfs)

    p17 = plot(rstar_s[4][1:50], -zc, lw=ls, lc=bc[4], label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    plot!(rstar_s[8][1:50], -zc, lw=ls, lc=bc[8], label="")
    plot!(rstar_s[12][1:50], -zc, lw=ls, lc=bc[12], label="")
    plot!(Ds[1:50, 4], -zc, lw=ls, lc=dc[4], linestyle=:dot, alpha=ab, label=" D4", legendfontsize=lfs)

    p18 = plot(bs[1:50, 4], -zc, lw=ls, lc=bc[4], label=" B4", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, legendfontsize=lfs)
    plot!(bs[1:50, 8], -zc, lw=ls, lc=bc[8], label=" B8", alpha=ab, legendfontsize=lfs)
    plot!(bs[1:50, 12], -zc, lw=ls, lc=bc[12], label=" B12", alpha=ab, legendfontsize=lfs)

    p19 = plot(rstar_s[5][1:50], -zc, lw=ls, lc=bc[5], label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    plot!(rstar_s[9][1:50], -zc, lw=ls, lc=bc[9], alpha=ab, label="")
    # plot!(rstar_s[13][1:50], -zc, lw=ls, lc=bc[13], label="")
    plot!(Ds[1:50, 5], -zc, lw=ls, lc=dc[5], linestyle=:dot, alpha=ab, label=" D5", legendfontsize=lfs)

    p20 = plot(bs[1:50, 5], -zc, lw=ls, lc=bc[5], label=" B5", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, legendfontsize=lfs)
    plot!(bs[1:50, 9], -zc, lw=ls, lc=bc[9], label=" B9", alpha=ab, legendfontsize=lfs)
    plot!(bs[1:50, 13], -zc, lw=ls, lc=bc[13], linestyle=:dash, label=" B13", alpha=ab_ext, legendfontsize=lfs)
    
    winter = plot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
    fg_legend = :transparent,
    layout = (1,10),
    size=(1600,900),
    plot_xaxis="Depth (m)"
    # xlabel = "R*",
    # plot_title="Winter", 
    # plot_titlefontsize = 20,
    # titlefontsize=tfs, titlelocation=:center, 
    )

    summer = plot(p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,
    fg_legend = :transparent,
    layout = (1,10),
    size=(1600,900),
    plot_xaxis="Depth (m)"
    # xlabel = "R*",
    # plot_title="Summer", 
    # plot_titlefontsize = 20,
    # titlefontsize=tfs, titlelocation=:center, 
    )
    
    combined = plot(winter, summer, 
        fg_legend = :transparent,
        layout = (2,1),
        size=(1600,900),
        # plot_xaxis="Depth (m)"
        # xlabel = "R*",
        # plot_title="Winter", 
        # titlefontsize=tfs, titlelocation=:center, 
    )

    savefig(combined,"/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rsB/rsB_$(fsaven).png")
    
    return combined

end

function plot_rstar_7B4D(fsaven, rstar_w, rstar_s, Dw, Ds, Bw, Bs)
    # B1 eats D1 (POM)  |  B2 & B7 eat D2  |  B3 & B6 eat D3  |  B4 & B5 eat D4

    H = 500
    dz = 10
    zc = [dz/2:dz:(H-dz/2)]
    lg=:bottomright
    tfs=22
    tl=:right

    #--------------------------------------------------------------------------------------------
    # Winter
    #--------------------------------------------------------------------------------------------
    p1 = plot(rstar_w[1][1:50], -zc, lw=ls, lc=bc[1], label="", ylabel="Depth (m)", xrotation=45, xlabel="", 
    border=:box, legend=lg, xscale=:log10, title="D1", title_loc=tl, titlefontsize=tfs)
    plot!(Dw[1:50, 1], -zc, lw=ls, lc=dc[1], linestyle=:dot, label=" D1", alpha=ab, legendfontsize=lfs)

    p2 = plot(Bw[1:50, 1], -zc, lw=ls, lc=bc[1], label=" B1", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xscale=:log10, legendfontsize=lfs)

    p3 = plot(rstar_w[3][1:50], -zc, lw=ls, lc=bc[3], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), title="D3", title_loc=tl, titlefontsize=tfs)
    # plot!(rstar_w[6][1:50], -zc, lw=ls, lc=bc[6],  label="")
    plot!(Dw[1:50, 3], -zc, lw=ls, lc=dc[3], linestyle=:dot, alpha=ab, label=" D3", legendfontsize=lfs)

    p4 = plot(Bw[1:50, 3], -zc, lw=ls, lc=bc[3], label=" B3", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, legendfontsize=lfs)
    plot!(Bw[1:50, 6], -zc, lw=ls, lc=bc[6],  label=" B6", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)

    p5 = plot(rstar_w[5][1:50], -zc, lw=ls, lc=bc[5], label="", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), title="D4", title_loc=tl, titlefontsize=tfs)
    # plot!(rstar_w[4][1:50], -zc, lw=ls, lc=bc[4], label="")
    plot!(Dw[1:50, 4], -zc, lw=ls, lc=dc[4], linestyle=:dot, alpha=ab, label=" D4", legendfontsize=lfs)

    p6 = plot(Bw[1:50, 4], -zc, lw=ls, lc=bc[4], label=" B4", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
    plot!(Bw[1:50, 5], -zc, lw=ls, lc=bc[5], label=" B5", alpha=ab, legendfontsize=lfs)

    #--------------------------------------------------------------------------------------------
    # Summer
    #--------------------------------------------------------------------------------------------
    p7 = plot(rstar_s[1][1:50], -zc, lw=ls, lc=bc[1], label="", ylabel="Depth (m)", 
    xrotation=45, xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10)
    plot!(Ds[1:50, 1], -zc, lw=ls, lc=dc[1], linestyle=:dot,label=" D1", alpha=ab, legendfontsize=lfs)

    p8 = plot(Bs[1:50, 1], -zc, lw=ls, lc=bc[1], label=" B1", xrotation=45, xguidefontsize=12, 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, xscale=:log10, xlabel=L"Biomass", legendfontsize=lfs)

    p9 = plot(rstar_s[6][1:50], -zc, lw=ls, lc=bc[6], label="", xrotation=45, xguidefontsize=12, 
    xlabel=L"R*", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), alpha=ab)
    # plot!(rstar_s[3][1:50], -zc, lw=ls, lc=bc[3], label="")
    plot!(Ds[1:50, 3], -zc, lw=ls, lc=dc[3], linestyle=:dot, alpha=ab, label=" D3", legendfontsize=lfs)

    p10 = plot(Bs[1:50, 3], -zc, lw=ls, lc=bc[3], label=" B3", xrotation=45, xguidefontsize=12, 
    xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
    plot!(Bs[1:50, 6], -zc, lw=ls, lc=bc[6], label=" B6", alpha=ab, legendfontsize=lfs)

    p11 = plot(rstar_s[4][1:50], -zc, lw=ls, lc=bc[4], label="", xrotation=45, xguidefontsize=12, xlabel=L"R*", 
    border=:box, legend=lg, xscale=:log10, yformatter=Returns(""))
    # plot!(rstar_s[5][1:50], -zc, lw=ls, lc=bc[5], label="")
    plot!(Ds[1:50, 4], -zc, lw=ls, lc=dc[4], linestyle=:dot, alpha=0.6, label=" D4", legendfontsize=lfs)

    p12 = plot(Bs[1:50, 4], -zc, lw=ls, lc=bc[4], label=" B4", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
    border=:box, legend=lg, yformatter=Returns(""), legendfontsize=lfs)
    plot!(Bs[1:50, 5], -zc, lw=ls, lc=bc[5], label=" B5", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
    
    
    combined = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
        fg_legend = :transparent,
        layout = (2,6),
        size=(900,700),
    )

    savefig(combined, "/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rsB/rsB_$(fsaven).png")

end

function plot_rstar_2B2D(fsaven, rstar_w, rstar_s, Dw, Ds, Bw, Bs)
    # B1 eats D1 (POM)  |  B2 eats D2 (DOM)  
    # zc = get_zc(ds, 500)
    H = 500
    dz = 10
    zc = [dz/2:dz:(H-dz/2)]
    lg=:bottomright
    tfs=12

    #--------------------------------------------------------------------------------------------
    # Winter
    #--------------------------------------------------------------------------------------------
    p1 = plot(rstar_w[1][1:50], -zc, lw=ls, lc=bc[1], ylabel="Depth (m)", label="", xrotation=45, 
    xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, alpha=ab, legendfontsize=lfs)
    plot!(Dw[1:50, 1], -zc, lw=ls, lc=dc[1], linestyle=:dot, label=" D1", alpha=ab, legendfontsize=lfs)

    p2 = plot(Bw[1:50, 1], -zc, lw=ls, lc=bc[1], label=" B1", xrotation=45, xguidefontsize=12, xlabel="", 
    border=:box, legend=lg, yformatter=Returns(""), alpha=ab, legendfontsize=lfs)

    p3 = plot(rstar_w[2][1:50], -zc, lw=ls, lc=bc[2], label="",xrotation=45, 
    xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), alpha=ab)
    plot!(Dw[1:50, 2], -zc, lw=ls, lc=dc[2], linestyle=:dot, alpha=ab, label=" D2", legendfontsize=lfs)

    p4 = plot(Bw[1:50, 2], -zc, lw=ls, lc=bc[2], ylabel="", label=" B2",xrotation=45, 
    xguidefontsize=12, xlabel="", border=:box, legend=lg, yformatter=Returns(""), alpha=ab, legendfontsize=lfs)


    #--------------------------------------------------------------------------------------------
    # Summer
    #--------------------------------------------------------------------------------------------
    p5 = plot(rstar_s[1][1:50], -zc, lw=ls, lc=bc[1], ylabel="Depth (m)", label="", xrotation=45, 
    xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, alpha=ab)
    plot!(Ds[1:50, 1], -zc, lw=ls, lc=dc[1], linestyle=:dot, alpha=ab, label=" D1", legendfontsize=lfs)

    p6 = plot(Bs[1:50, 1], -zc, lw=ls, lc=bc[1], label=" B1", xrotation=45, xguidefontsize=12, 
    border=:box, legend=lg, yformatter=Returns(""), xlabel=L"Biomass", alpha=ab, legendfontsize=lfs)

    p7 = plot(rstar_s[2][1:50], -zc, lw=ls, lc=bc[2], label="",xrotation=45, 
    xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, yformatter=Returns(""), alpha=ab)
    plot!(Ds[1:50, 2], -zc, lw=ls, lc=dc[2], linestyle=:dot, alpha=ab, label=" D2", legendfontsize=lfs)

    p8 = plot(Bs[1:50, 2], -zc, lw=ls, lc=bc[2], ylabel="", label=" B2",xrotation=45, 
    xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, yformatter=Returns(""), alpha=ab, legendfontsize=lfs)

    combined = plot(p1, p2, p3, p4, p5, p6, p7, p8,
        fg_legend = :transparent,
        fontfamily="Computer Modern",
        layout = (2,4),
        size=(600,500),
    )

    savefig(combined, "/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rsB/rsB_$(fsaven).png")

end

function plot_rstar_P(fsaven, rstar_w, rstar_s, Pw, Ps, Nw, Ns, np)

    H = 200
    dz = 10
    zc = [dz/2:dz:(H-dz/2)]
    lg=:bottomright
    tfs=12
    lfs=9


    if np == 1
        p1 = plot(rstar_w[1][1:20], -zc, lw=ls, lc=pc[1], label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, title="", alpha=ab)
        plot!(Nw[1:20, 1], -zc, lw=ls, lc=nc[1], linestyle=:dot, label="", legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=ls, lc=pc[1], label="", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), alpha=ab)

        p3 = plot(rstar_s[1][1:20], -zc, lw=ls, lc=pc[1], label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, title="", alpha=ab)
        plot!(Ns[1:20, 1], -zc, lw=ls, lc=nc[1], linestyle=:dot,label=" N", legend=:bottom, legendfontsize=lfs, alpha=ab)

        p4 = plot(Ps[1:20, 1], -zc, lw=ls, lc=pc[1], label=" P1", ylabel="", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
        border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), legendfontsize=lfs, alpha=ab)
    end
    
    if np == 2
        p1 = plot(rstar_w[1][1:20], -zc, lw=ls, lc=pc[1], label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(rstar_w[2][1:20], -zc, lw=ls, lc=pc[2],  label="")
        plot!(Nw[1:20, 1], -zc, lw=ls, lc=nc[1], linestyle=:dot, label="", alpha=ab, legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=ls, lc=pc[1], label="", ylabel="", xrotation=45, alpha=ab, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""))
        plot!(Pw[1:20, 2], -zc, lw=ls, lc=pc[2],  label="", alpha=ab)

        p3 = plot(rstar_s[1][1:20], -zc, lw=ls, lc=pc[1], label="", ylabel="Depth (m)", xrotation=45, alpha=ab, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        plot!(rstar_s[2][1:20], -zc, lw=ls, lc=pc[2], alpha=ab,  label="")
        plot!(Ns[1:20, 1], -zc, lw=ls, lc=nc[1], linestyle=:dot,label=" N", alpha=ab, legend=:bottom, legendfontsize=lfs)

        p4 = plot(Ps[1:20, 1], -zc, lw=ls, lc=pc[1], label=" P1", ylabel="", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
        border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), alpha=ab, legendfontsize=lfs)
        plot!(Ps[1:20, 2], -zc, lw=ls, lc=pc[2],  label=" P2", alpha=ab, legendfontsize=lfs)
    end 

    if np == 4
        p1 = plot(rstar_w[2][1:20], -zc, lw=ls, lc=pc[2], label="", ylabel="Depth (m)", xrotation=45, alpha=ab, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10)
        # plot!(rstar_w[1][1:20], -zc, lw=ls, lc=pc[1], alpha=ab,  label="")
        # plot!(rstar_w[3][1:20], -zc, lw=ls, lc=pc[3], alpha=ab,  label="")
        plot!(rstar_w[4][1:20], -zc, lw=ls, lc=pc[4], alpha=ab,  label="")
        plot!(Nw[1:20, 1], -zc, lw=ls, lc=nc[1], linestyle=:dot, alpha=ab, label="", legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=ls, lc=pc[1], label="", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext)
        plot!(Pw[1:20, 2], -zc, lw=ls, lc=pc[2], alpha=ab,  label="")
        plot!(Pw[1:20, 3], -zc, lw=ls, lc=pc[3], label="", linestyle=:dash, alpha=ab_ext)
        plot!(Pw[1:20, 4], -zc, lw=ls, lc=pc[4], alpha=ab,  label="")

        p3 = plot(rstar_s[1][1:20], -zc, lw=ls, lc=pc[1], label="", ylabel="Depth (m)", xrotation=45, alpha=ab, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, title="", titlefontsize=tfs)
        # plot!(rstar_s[2][1:20], -zc, lw=ls, lc=pc[2],  label="")
        plot!(rstar_s[3][1:20], -zc, lw=ls, lc=pc[3], alpha=ab, label="")
        # plot!(rstar_s[4][1:20], -zc, lw=ls, lc=pc[4],  label="")
        plot!(Ns[1:20, 1], -zc, lw=ls, lc=nc[1], linestyle=:dot,label=" N", legend=:bottom, alpha=ab, legendfontsize=lfs)

        p4 = plot(Ps[1:20, 1], -zc, lw=ls, lc=pc[1], label=" P1", ylabel="", xrotation=45, alpha=ab, 
        xguidefontsize=12, xlabel=L"Biomass", border=:box, legend=lg, title="", titlefontsize=tfs, yformatter=Returns(""), legendfontsize=lfs)
        plot!(Ps[1:20, 2], -zc, lw=ls, lc=pc[2],  label=" P2", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
        plot!(Ps[1:20, 3], -zc, lw=ls, lc=pc[3],  label=" P3", alpha=ab, legendfontsize=lfs)
        plot!(Ps[1:20, 4], -zc, lw=ls, lc=pc[4],  label=" P4", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs) 
    end

    if np == 8
        p1 = plot(rstar_w[2][1:20], -zc, lw=ls, lc=pc[2], label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, xscale=:log10, alpha=ab)
        # plot!(rstar_w[1][1:20], -zc, lw=ls, lc=pc[1],  label="", alpha=ab)
        # plot!(rstar_w[3][1:20], -zc, lw=ls, lc=pc[3],  label="", alpha=ab)
        # plot!(rstar_w[4][1:20], -zc, lw=ls, lc=pc[4],  label="", alpha=ab)
        # plot!(rstar_w[5][1:20], -zc, lw=ls, lc=pc[5],  label="", alpha=ab)
        plot!(rstar_w[6][1:20], -zc, lw=ls, lc=pc[6],  label="", alpha=ab)
        # plot!(rstar_w[7][1:20], -zc, lw=ls, lc=pc[7],  label="", alpha=ab)
        # plot!(rstar_w[8][1:20], -zc, lw=ls, lc=pc[8],  label="", alpha=ab)
        plot!(Nw[1:20, 1], -zc, lw=ls, lc=nc[1], linestyle=:dot, label="", alpha=ab, legend=:bottom)

        p2 = plot(Pw[1:20, 1], -zc, lw=ls, lc=pc[1], label="", ylabel="", xrotation=45, 
        xguidefontsize=12, xlabel="", border=:box, legend=lg, linestyle=:dash, alpha=ab_ext,  yformatter=Returns(""), xlims=(-0.01, 0.15))
        plot!(Pw[1:20, 2], -zc, lw=ls, lc=pc[2], alpha=ab, label="")
        plot!(Pw[1:20, 3], -zc, lw=ls, lc=pc[3], linestyle=:dash, alpha=ab_ext, label="")
        plot!(Pw[1:20, 4], -zc, lw=ls, lc=pc[4], linestyle=:dash, alpha=ab_ext, label="")
        plot!(Pw[1:20, 5], -zc, lw=ls, lc=pc[5], linestyle=:dash, alpha=ab_ext, label="")
        plot!(Pw[1:20, 6], -zc, lw=ls, lc=pc[6], alpha=ab,  label="")
        plot!(Pw[1:20, 7], -zc, lw=ls, lc=pc[7], linestyle=:dash, alpha=ab_ext, label="")
        plot!(Pw[1:20, 8], -zc, lw=ls, lc=pc[8], linestyle=:dash, alpha=ab_ext, label="")

        p3 = plot(rstar_s[4][1:20], -zc, lw=ls, lc=pc[4], label="", ylabel="Depth (m)", xrotation=45, 
        xguidefontsize=12, xlabel=L"R*", border=:box, legend=lg, xscale=:log10, alpha=ab)
        # plot!(rstar_s[1][1:20], -zc, lw=ls, lc=pc[1],  label="", alpha=ab)
        # plot!(rstar_s[2][1:20], -zc, lw=ls, lc=pc[2],  label="", alpha=ab)
        # plot!(rstar_s[3][1:20], -zc, lw=ls, lc=pc[3],  label="", alpha=ab)
        # plot!(rstar_s[4][1:20], -zc, lw=ls, lc=pc[4],  label="", alpha=ab)
        plot!(rstar_s[5][1:20], -zc, lw=ls, lc=pc[5],  label="", alpha=ab)
        # plot!(rstar_s[6][1:20], -zc, lw=ls, lc=pc[6],  label="", alpha=ab)
        # plot!(rstar_s[7][1:20], -zc, lw=ls, lc=pc[7],  label="", alpha=ab)
        # plot!(rstar_s[8][1:20], -zc, lw=ls, lc=pc[8],  label="", alpha=ab)
        plot!(Ns[1:20, 1], -zc, lw=ls, lc=nc[1], linestyle=:dot,label=" N", legend=:bottom, alpha=ab, legendfontsize=lfs)

        p4 = plot(Ps[1:20, 1], -zc, lw=ls, lc=pc[1], label=" P1", ylabel="", xrotation=45, xguidefontsize=12, xlabel=L"Biomass", 
        border=:box, legend=lg, yformatter=Returns(""), linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
        plot!(Ps[1:20, 2], -zc, lw=ls, lc=pc[2],  label=" P2", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
        plot!(Ps[1:20, 3], -zc, lw=ls, lc=pc[3],  label=" P3", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
        plot!(Ps[1:20, 4], -zc, lw=ls, lc=pc[4],  label=" P4", legendfontsize=lfs, alpha=ab ) 
        plot!(Ps[1:20, 5], -zc, lw=ls, lc=pc[5],  label=" P5", legendfontsize=lfs, alpha=ab)
        plot!(Ps[1:20, 6], -zc, lw=ls, lc=pc[6],  label=" P6", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs)
        plot!(Ps[1:20, 7], -zc, lw=ls, lc=pc[7],  label=" P7", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs) 
        plot!(Ps[1:20, 8], -zc, lw=ls, lc=pc[8],  label=" P8", linestyle=:dash, alpha=ab_ext, legendfontsize=lfs) 

    end
    

    combined = plot(p1, p2, p3, p4,
    fg_legend = :transparent,
    layout = (2,2),
    size=(450,600),
    plot_title=L"R* P_N"
    )

    savefig(combined,"/home/lee/Dropbox/Development/NPZBD_1D/results/plots/rstar/rsP/rsP_$(fsaven).png")

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
        RS[i] = ifelse(RS[i] > 15, NaN, RS[i])
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
# - ws[1] = 6.0 
# win_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230915_00:03_8P6Z13B5D_ep.nc"
# sum_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Su100y_230915_17:36_8P6Z13B5D_ep.nc"

# - ws[1] = 10.0 
# win_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230915_22:19_8P6Z13B5D_ep.nc"
# sum_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Su100y_230916_11:25_8P6Z13B5D_ep.nc"

# winter = NCDataset(win_p)
# summer = NCDataset(sum_p)
# rgx = match(r"(/[A-Z][a-z])(?=\d)", win_p)
# wp = replace(win_p, rgx[1] => "/WiSu")
# fsaven = replace(wp, "_ep.nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/" => "")

# CMd = display(winter["CM"][:])
# GrMd = display(winter["GrM"][:])

# # Get endpoints and set all biomass values below cut-off (10^-6) equal to cut-off
# Nw, Pw, Zw, Bw, Dw = get_endpoints(["n", "p", "z", "b", "d"], winter)
# Ns, Ps, Zs, Bs, Ds = get_endpoints(["n", "p", "z", "b", "d"], summer)

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

# plot_rstar_13B5D(fsaven, rstarB_w, rstarB_s, Dw, Ds, Bwe, Bse)
# plot_rstar_P(fsaven, rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, 8)

# # -----------------------------------------------------------------------
# # 1N 4P 3Z 7B 4D
# # -----------------------------------------------------------------------
# win_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230827_13:45_4P3Z7B4D_ep.nc"
# sum_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Su100y_230827_17:10_4P3Z7B4D_ep.nc"

# rgx = match(r"(/[A-Z][a-z])(?=\d)", win_p)
# wp = replace(win_p, rgx[1] => "/WiSu")
# fsaven = replace(wp, "_ep.nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/" => "")

# winter = NCDataset(win_p)
# summer = NCDataset(sum_p)
# CMd = display(winter["CM"][:])
# GrMd = display(winter["GrM"][:])

# # Get endpoints and set all biomass values below cut-off (10^-6) equal to cut-off
# Nw, Pw, Zw, Bw, Dw = get_endpoints(["n", "p", "z", "b", "d"], winter)
# Ns, Ps, Zs, Bs, Ds = get_endpoints(["n", "p", "z", "b", "d"], summer)

# Bww = cut_off(Bw, 7)
# Bss = cut_off(Bs, 7)
# Pww = cut_off(Pw, 4)
# Pss = cut_off(Ps, 4)

# # Calculate Rstar B
# rstarB_ijw = get_rstar_B(Bww, Zw, winter, 7, 3, "Win")
# rstarB_ijs = get_rstar_B(Bss, Zs, summer, 7, 3, "Sum")
# rstarP_w = get_rstar_P(Pww, Zw, winter, 4, 3, "Win")
# rstarP_s = get_rstar_P(Pss, Zs, summer, 4, 3, "Sum")

# Bwe = extinct(Bww, 7)
# Bse = extinct(Bss, 7)
# Pwe = extinct(Pww, 4)
# Pse = extinct(Pss, 4)

# plot_rstar_7B4D(fsaven, rstarB_ijw, rstarB_ijs, Dw, Ds, Bwe, Bse)
# plot_rstar_P(fsaven, rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, 4)
# println("done")


# # -----------------------------------------------------------------------
# # 1N 2P 2Z 2B 2D
# # -----------------------------------------------------------------------
# win_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230905_20:05_2P2Z2B2D_ep.nc"
# sum_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Su100y_230905_23:01_2P2Z2B2D_ep.nc"

# rgx = match(r"(/[A-Z][a-z])(?=\d)", win_p)
# wp = replace(win_p, rgx[1] => "/WiSu")
# fsaven = replace(wp, "_ep.nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/" => "")

# winter = NCDataset(win_p)
# summer = NCDataset(sum_p)
# CMd = display(winter["CM"][:])
# GrMd = display(winter["GrM"][:])

# # Get endpoints 
# Nw, Pw, Zw, Bw, Dw = get_endpoints(["n", "p", "z", "b", "d"], winter)
# Ns, Ps, Zs, Bs, Ds = get_endpoints(["n", "p", "z", "b", "d"], summer)

# # set all biomass values below cut-off (10^-6) equal to cut-off
# Bww = cut_off(Bw, 2)
# Bss = cut_off(Bs, 2)
# Pww = cut_off(Pw, 2)
# Pss = cut_off(Ps, 2)

# # Calculate Rstar 
# rstarB_ijw = get_rstar_B(Bww, Zw, winter, 2, 2, "Win")
# rstarB_ijs = get_rstar_B(Bss, Zs, summer, 2, 2, "Sum")
# rstarP_w = get_rstar_P(Pww, Zw, winter, 2, 2, "Win")
# rstarP_s = get_rstar_P(Pss, Zs, summer, 2, 2, "Sum")

# # set low biomass to zero
# Bwe = extinct(Bww, 2)
# Bse = extinct(Bss, 2)
# Pwe = extinct(Pww, 2)
# Pse = extinct(Pss, 2)


# # Plot
# plot_rstar_2B2D(fsaven, rstarB_ijw, rstarB_ijs, Dw, Ds, Bwe, Bse)
# plot_rstar_P(fsaven, rstarP_w, rstarP_s, Pwe, Pse, Nw, Ns, 2)
# println("done")


# # -----------------------------------------------------------------------
# # 1N 1P 3Z 2B 2D
# # -----------------------------------------------------------------------
win_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Wi100y_230906_13:51_1P3Z2B2D_ep.nc"
sum_p = "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/Su100y_230906_14:13_1P3Z2B2D_ep.nc"

rgx = match(r"(/[A-Z][a-z])(?=\d)", win_p)
wp = replace(win_p, rgx[1] => "/WiSu")
fsaven = replace(wp, "_ep.nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/endpoints/" => "")

winter = NCDataset(win_p)
summer = NCDataset(sum_p)
CMd = display(winter["CM"][:])
GrMd = display(winter["GrM"][:])

# Get endpoints
Nw, Pw, Zw, Bw, Dw = get_endpoints(["n", "p", "z", "b", "d"], winter)
Ns, Ps, Zs, Bs, Ds = get_endpoints(["n", "p", "z", "b", "d"], summer)

# set all biomass values below cut-off (10^-6) equal to cut-off
Bww = cut_off(Bw, 2)
Bss = cut_off(Bs, 2)
Pww = cut_off(Pw, 1)
Pss = cut_off(Ps, 1)

# Calculate Rstar 
rstarBw = get_rstar_B(Bww, Zw, winter, 2, 3, "Win")
rstarBs = get_rstar_B(Bss, Zs, summer, 2, 3, "Sum")
rstarPw = get_rstar_P(Pww, Zw, winter, 1, 3, "Win")
rstarPs = get_rstar_P(Pss, Zs, summer, 1, 3, "Sum")

Bwe = extinct(Bww, 2)
Bse = extinct(Bss, 2)
Pwe = extinct(Pww, 1)
Pse = extinct(Pss, 1)

# Plot
plot_rstar_2B2D(fsaven, rstarBw, rstarBs, Dw, Ds, Bwe, Bse)
plot_rstar_P(fsaven, rstarPw, rstarPs, Pwe, Pse, Nw, Ns, 1)
println("done")



