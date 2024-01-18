using KernelDensity, Plots, LaTeXStrings


function ridgeline_plot(data, plotxlab, colname, scaling::Bool=false)

    #----------------------------------------------------------
    # plotxlab = L"\degree C"
    plottitle=""
    plotylab = ""
    ylabels = ["Winter", "Spring", "Summer", "Fall"]
    ridgecolors = ["blue", "seagreen4", "red", "yellow3"]
    spacer = 0.5
    riser = 0.001
    xls = 10
    yls = 10
    tls = 15
    hlinecolor = "black"
    hlw = 0.5
    halpha = 0.6
    ridgealpha = 0.9
    ridgeoutline = "black"
    ridgelw = 0.5
    showplot::Bool = false
    #----------------------------------------------------------

    # scaling (to show relavtive densities)
    if scaling == true
        scaling = length.(data) ./ maximum(length.(data))
    else
        scaling = ones(length(data))
    end

    # make density plot line data for plotting
    dense, xs = Any[], Any[]
    for i in 1:size(data, 1)
        temp = kde(data[i])
        push!(dense, [(temp.density .* scaling[i]) .+ ((((i - 1) * spacer) .+ riser))])
        push!(xs, [temp.x])
    end

    # find ideal x and y lims
    xlimits, ylimits = Float64[], Float64[]
    for i in 1:size(data, 1)
        push!(xlimits, minimum(xs[i][1]))
        push!(xlimits, maximum(xs[i][1]))
        push!(ylimits, minimum(dense[i][1]))
        push!(ylimits, maximum(dense[i][1]))
    end

    colname == "CTDFLUOR" ? xl=(minimum(xlimits), 5) : xl=[minimum(xlimits), maximum(xlimits)]

    # plot fig
    fig = plot(yticks = (collect((size(dense,1) - 1):-1:0) .* (spacer), reverse(ylabels)))
        plot!(grid=false, border=:box, title=plottitle, xlab=plotxlab, xrotation=45, ylab = plotylab,
        xtickfontsize=xls, ytickfontsize=yls, titlefontsize=tls, 
        ylim= [minimum(ylimits), maximum(ylimits)], xlim=xl)
        hline!([(collect((size(dense,1) - 1):-1:0)) .* (spacer)], color = hlinecolor, lw = hlw, label = "", alpha = halpha)


    #plotting each curve
    for i in size(dense, 1):-1:1

        #add plots
        if showplot == true
            display(plot!(xs[i], dense[i], fillrange = ((i - 1) * spacer), label = "", fillalpha = ridgealpha,
                    linecolor = ridgeoutline, lw = ridgelw, fillcolor = ridgecolors[i]))
        else
            plot!(xs[i], dense[i], fillrange = ((i - 1) * spacer), label = "", fillalpha = ridgealpha,
                    linecolor = ridgeoutline, lw = ridgelw, fillcolor = ridgecolors[i]);
        end

    end

    return(fig)

end


function plot_hist2d(wi, sp, su, fa, overall_mean, xlab, colname)

    lfs=8
    tfs=10
    lwmean=3
    lcmean="cyan"
    alp=0.8
    bin=30

    if colname == "CTDFLUOR"
        histlg=:topright 
    elseif colname == "CTDOXY" || colname == "CTDATTN"
        histlg=:bottomright 
    else 
        histlg=:topleft
    end

    colname == "CTDFLUOR" ? xl=(0, 6.5) : xl=()

    p1 = histogram2d(wi[!,colname], -wi[!,"CTDDEPTH"], bins=bin, grid=false, ylabel="", xlabel="", 
                    xrotation=45, colorbar=false, title = "Winter", border=:box, alpha=alp, titlefontsize=tfs)
                    plot!(overall_mean[!, "all_mean"], -overall_mean[!,"CTDDEPTH"], lw=lwmean, lc=lcmean, label=false)
    p2 = histogram2d(sp[!,colname], -sp[!,"CTDDEPTH"], bins=bin, grid=false, ylabel="Depth (m)", xlabel="", 
                    xrotation=45, colorbar=false, ymirror = true, title = "Spring", border=:box, alpha=alp, titlefontsize=tfs)
                    plot!(overall_mean[!, "all_mean"], -overall_mean[!,"CTDDEPTH"], lw=lwmean, lc=lcmean, label=false)
    p3 = histogram2d(su[!,colname], -su[!,"CTDDEPTH"], bins=bin, grid=false, ylabel="", xlabel=xlab, 
                    xrotation=45, colorbar=false, title = "Summer", border=:box, alpha=alp, titlefontsize=tfs)
                    plot!(overall_mean[!, "all_mean"], -overall_mean[!,"CTDDEPTH"], lw=lwmean, lc=lcmean, label=false)
    p4 = histogram2d(fa[!,colname], -fa[!,"CTDDEPTH"], bins=bin, grid=false, ylabel="Depth (m)", xlabel=xlab, 
                    xrotation=45, colorbar=false, ymirror = true, title = "Fall", border=:box, alpha=alp, titlefontsize=tfs)
                    plot!(overall_mean[!, "all_mean"], -overall_mean[!,"CTDDEPTH"], lw=lwmean, lc=lcmean, label=" Overall\n Mean", 
                    legend=histlg, legendfontsize=lfs)

    return p1, p2, p3, p4

end


function plot_line_box_violin(wi, sp, su, fa, wi_means, sp_means, su_means, fa_means, xlab, colname, maxdepth)

    lfs=8
    maxdepth == 0 ? meanlw=8 : meanlw=12
    meanalph=0.6
    maxdepth == 0 ? lgloc=:topleft : lgloc=:bottomright
    colname == "CTDFLUOR" ? xl=(0, 6.5) : xl=()

    x = [["Winter"], ["Spring"], ["Summer"], ["Fall"]]
    y = [[wi[!,colname]], [sp[!,colname]], [su[!,colname]], [fa[!,colname]]]

    p6 = plot(wi_means[!,colname], -wi_means[!,"CTDDEPTH"], lw=meanlw, lc=:blue, label=" Winter", legendfontsize=lfs, legend=lgloc, 
              foreground_color_legend = nothing, legend_title=" Means", legend_title_font_halign=:right, legend_background_color=:transparent,
              ylabel="Depth (m)", xlabel=xlab, grid=false, border=:box, xrotation=45, alpha=meanalph)
         plot!(sp_means[!,colname], -sp_means[!,"CTDDEPTH"], lw=meanlw, lc=:seagreen4, label=" Spring", legendfontsize=lfs, legend=lgloc, alpha=meanalph)
         plot!(su_means[!,colname], -su_means[!,"CTDDEPTH"],lw=meanlw, lc=:red, label=" Summer", legendfontsize=lfs, legend=lgloc, alpha=meanalph)
         plot!(fa_means[!,colname], -fa_means[!,"CTDDEPTH"],lw=meanlw, lc=:yellow3, label=" Fall", legendfontsize=lfs, legend=lgloc, alpha=meanalph)

    if maxdepth == 0

        boxplot!(x[1], y[1] ,grid=false, label=false, line=(2, :black), fill=(0.3, :blue), inset=(1, bbox(0.03, 0.15, 0.45, 0.7, :top, :right)), 
                subplot=2, border=:box, bg_color_inside=:snow, permute=(:x, :y), outliers=false)
        violin!(x[1], y[1], grid=false, label=false, line=(0.5, :black), fill=(0.5, :grey), subplot=2, permute=(:x, :y), xlim=xl)
        boxplot!(x[2], y[2], grid=false, label=false, line=(2, :black), fill=(0.3, :seagreen4), subplot=2, permute=(:x, :y), outliers=false)
        violin!(x[2], y[2], grid=false, label=false, line=(0.5, :black), fill=(0.5, :grey), subplot=2, permute=(:x, :y), xlim=xl)
        boxplot!(x[3], y[3], grid=false, label=false, line=(2, :black), fill=(0.3, :red), subplot=2, permute=(:x, :y), outliers=false)
        violin!(x[3], y[3], grid=false, label=false, line=(0.5, :black), fill=(0.5, :grey), subplot=2, permute=(:x, :y), xlim=xl)
        boxplot!(x[4], y[4], grid=false, label=false, line=(2, :black), fill=(0.3, :yellow3), subplot=2, permute=(:x, :y), outliers=false)
        violin!(x[4], y[4], grid=false, label=false, line=(0.5, :black), fill=(0.5, :grey), subplot=2, permute=(:x, :y), xlim=xl)

        p7 = nothing

    else

        xls, yls = 10, 10  
        p7 =  boxplot(x[1], y[1] ,grid=false, label=false, line=(2, :black), fill=(0.3, :blue), border=:box, permute=(:x, :y), 
            outliers=false, yrotation=45, xmirror = true, xtickfontsize=xls, ytickfontsize=yls, ylabel=xlab)
            violin!(x[1], y[1], grid=false, label=false, line=(0.5, :black), fill=(0.5, :grey), permute=(:x, :y))
            boxplot!(x[2], y[2], grid=false, label=false, line=(2, :black), fill=(0.3, :seagreen4), permute=(:x, :y), outliers=false)
            violin!(x[2], y[2], grid=false, label=false, line=(0.5, :black), fill=(0.5, :grey), permute=(:x, :y))
            boxplot!(x[3], y[3], grid=false, label=false, line=(2, :black), fill=(0.3, :red), permute=(:x, :y), outliers=false)
            violin!(x[3], y[3], grid=false, label=false, line=(0.5, :black), fill=(0.5, :grey), permute=(:x, :y))
            boxplot!(x[4], y[4], grid=false, label=false, line=(2, :black), fill=(0.3, :yellow3), permute=(:x, :y), outliers=false)
            violin!(x[4], y[4], grid=false, label=false, line=(0.5, :black), fill=(0.5, :grey), permute=(:x, :y))

    end


    return p6, p7

end



# @userplot MarginalHist

# @recipe function f(h::MarginalHist)
#     if length(h.args) != 2 || !(typeof(h.args[1]) <: AbstractVector) ||
#         !(typeof(h.args[2]) <: AbstractVector)
#         error("Marginal Histograms should be given two vectors.  Got: $(typeof(h.args))")
#     end
#     x, y = h.args

#     # set up the subplots
#     legend := false
#     link := :both
#     framestyle := [:none :axes :none]
#     grid := false
#     layout := @layout [tophist           _
#                        hist2d{0.9w,0.9h} _]

#     # main histogram2d
#     @series begin
#         seriestype := :histogram2d
#         subplot := 2
#         x, y
#     end

#     # these are common to both marginal histograms
#     fillcolor := :black
#     fillalpha := 0.5
#     linealpha := 0.3
#     seriestype := :histogram
#     # margin = -10Plots.px

#     # upper histogram
#     @series begin
#         subplot := 1
#         bottom_margin = -10Plots.px
#         x
#     end

#     # right histogram
#     # @series begin
#     #     orientation := :h
#     #     subplot := 3
#     #     y
#     # end
# end