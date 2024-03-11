using DataFrames, CSV, Dates, StatsPlots, Statistics, CategoricalArrays, Measures
include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")

function load_guilds_csv(path)

    guilds_csv = DataFrame(CSV.File(path))

    return guilds_csv

end

function load_monthly_csvs(path, mean=false)

    month_names = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    out = Any[]
    sem_out = Any[]

    if mean == false
        for m in month_names
            csv = DataFrame(CSV.File("$(path)/$(m).csv"))
            push!(out, csv)
        end
    else
        for m in month_names
            csv = DataFrame(CSV.File("$(path)/RA_FL/$(m)_mean_RA_FL.csv"))
            sem = DataFrame(CSV.File("$(path)/sem_FL/$(m)_mean_RA_FL_sem.csv"))
            push!(out, csv)
            push!(sem_out, sem)
        end
    end

    return out, sem_out

end 

function prepare_groupedbar_inputs(month_names, data_matrix,  guilds)
      
    # Redefine unique for `CategoricalArray` types to return a categorical array, rather than a regular vector/array. 
    @eval function Base.unique(ctg::CategoricalArray) # can be run in REPL instead
        l = levels(ctg)
        newctg = CategoricalArray(l)
        levels!(newctg, l)
    end
  
    plot_months = repeat(month_names, inner = size(data_matrix)[1])
    plot_guilds = repeat(guilds, outer = size(data_matrix)[2])

    plot_mons = categorical(plot_months; levels = month_names)
    plot_glds = categorical(plot_guilds; levels = guilds)
  
    return plot_mons, data_matrix, plot_glds
  
end

function create_data_matrices(monthly_means, depth_rownums, season_nums, clades)
    # depth_rownums -> 1=5m, 2=DCM, 3=150m, 4=500m, 5=890m
    # returns matrix for each requested depth - cols are guild, rows are depth 

    clades = ["SAR11_clade", "Flavobacteriales", "SAR324_clade", "SAR202_clade"]
    # clades = ["SAR11_clade", "Flavobacteriales", "Ammonia_oxidizers", "Thioglobaceae"]
    
    depth_matrices = Any[]
    for i in eachindex(depth_rownums)
    
        month1 = monthly_means[season_nums[1]][depth_rownums[i], clades]
        month2 = monthly_means[season_nums[2]][depth_rownums[i], clades]
        month3 = monthly_means[season_nums[3]][depth_rownums[i], clades]
        season_df = DataFrame(month1) ; push!(season_df, month2); push!(season_df, month3)

        m = Matrix(season_df)
        mm = Matrix(transpose(m))
        push!(depth_matrices, mm)
    end

    return depth_matrices

end

function plot_grouped_bar_means(depth_matrices, depthstrings, season, month_names, bgcol)

    clade_labels = [" SAR11", " Flavo", " SAR324", " SAR202"];
    # clade_labels = [" SAR11", " Flavo", " Ammo", " Thio"];
    
    grouped_inputs = Any[]
    for i in eachindex(depth_matrices)
        mons, matrix, guilds = prepare_groupedbar_inputs(month_names, depth_matrices[i], clade_labels);
        push!(grouped_inputs, [mons, matrix, guilds])
    end

    yl=(0, 0.45)
    a=0.7
    lfs = 10
    xtfs = 12
    tfs=22

    plts_arr = Any[]
    if season=="Winter"

        for i in range(1, 5)
            if i == 1
                p_i = groupedbar(grouped_inputs[1][1], grouped_inputs[1][2], group=grouped_inputs[1][3],
                title = "$(depthstrings[1])\n",
                xlabel = "", ylabel = "Relative Abundance",
                legend=false, border=:box,
                grid=false, xrotation=45,
                ylims=yl, alpha=a,
                xtickfontsize=xtfs, 
                left_margin=5mm,
                titlefontsize=tfs,
                # background_inside=bgcol,
                )
            elseif i == 5
                p_i = groupedbar(grouped_inputs[5][1], grouped_inputs[5][2], group=grouped_inputs[5][3],
                title = "$(depthstrings[5])\n",
                xlabel = "", yformatter=Returns(""),
                legend_position = :top, border=:box,
                foreground_color_legend = nothing,
                background_color_legend = nothing,
                legendfontsize=lfs,
                grid=false, xrotation=45,  
                ylims=yl, alpha=a,
                xtickfontsize=xtfs, 
                titlefontsize=tfs,
                # background_inside=bgcol,
                )
            else
                p_i = groupedbar(grouped_inputs[i][1], grouped_inputs[i][2], group=grouped_inputs[i][3],
                title = "$(depthstrings[i])\n",
                xlabel = "", yformatter=Returns(""),
                legend=false, border=:box,
                grid=false, xrotation=45,
                ylims=yl, alpha=a,
                legendfontsize=lfs,
                xtickfontsize=xtfs, 
                titlefontsize=tfs,
                # background_inside=bgcol,
                )
            end
            push!(plts_arr, p_i)
        end

    else

        for i in range(1, 5)
            if i == 1
                p_i = groupedbar(grouped_inputs[1][1], grouped_inputs[1][2], group=grouped_inputs[1][3],
                xlabel = "", ylabel = "Relative Abundance",
                legend=false, border=:box,
                grid=false, xrotation=45,
                ylims=yl, alpha=a,
                xtickfontsize=xtfs, 
                left_margin=5mm,
                # background_inside=bgcol
                )
            else           
                p_i = groupedbar(grouped_inputs[i][1], grouped_inputs[i][2], group=grouped_inputs[i][3],
                xlabel = "", yformatter=Returns(""),
                legend=false, border=:box,
                grid=false, xrotation=45,
                ylims=yl, alpha=a,
                legendfontsize=lfs,
                xtickfontsize=xtfs, 
                titlefontsize=tfs,
                # background_inside=bgcol,
                )
            end
            push!(plts_arr, p_i)
        end
    end

    return plts_arr[1], plts_arr[2], plts_arr[3], plts_arr[4], plts_arr[5]

end


function plot_gbars(mpath)

    depths=["5m", "DCM", "150m", "500m", "890m"]
    season_names = ["Winter", "Spring", "Summer", "Fall"];
    month_names = [["Dec", "Jan", "Feb"], ["Mar", "Apr", "May"], ["Jun", "Jul", "Aug"], ["Sep", "Oct", "Nov"]];
    season_nums = [[12,1,2], [3,4,5], [9,10,11], [9,10,11]];
    depth_rownums = [1,2,3,4,5];
    bg_cols=["ghostwhite", "honeydew", "lavenderblush","linen"]

    monthly_means = load_monthly_csvs(mpath, 1);
    season_matrices = Any[]
    for s in eachindex(season_nums)
        matrices = create_data_matrices(monthly_means, depth_rownums, season_nums[s])
        push!(season_matrices, matrices)
    end

    season_plots = Any[]
    bg=1
    for i in eachindex(season_names)
        plts = plot_grouped_bar(season_matrices[i], depths, season_names[i], month_names[i], bg_cols[bg])
        push!(season_plots, plts...)
        bg+=1
    end

    fig = plot(season_plots..., layout=(4,5), size=(1000,1000))
    savefig(fig, "results/plots/SPOT/monthly_RA/monthly_RA_gbars_1.pdf")
    
end

#--------------------------------------------------------------------------------------
function filter_by_year(month_df, yr)

    data_by_year = month_df[in([yr]).(Dates.year.(month_df.Date)), :];

    return data_by_year

end ;

function filter_by_depth(month_df, dep, mean=false)

    if mean==false
        data_by_depth = filter(row -> row.Depth == dep, month_df)
    else
        data_by_depth = filter(row -> row.DepthBin == dep, month_df)
    end

    return data_by_depth

end ;

function get_year_depth_df(all_months_data, yr, dep, month_nums, mean=false)
       
    out = Any[]
    if mean==false
        for num in month_nums
            data = filter_by_depth(filter_by_year(all_months_data[num], yr), dep);
            push!(out, data)
        end
    else
        for num in month_nums
            data = filter_by_depth(all_months_data[num], dep, mean);
            push!(out, data)
        end
    end

    merged_df = vcat(out...)

    return merged_df
end ;


function create_data_matrices_months(months_dep_df, all_clade_colnames)

    mat = Matrix(months_dep_df[!, all_clade_colnames])

    return mat

end ;

function prep_gbar_bloom(month_names, guild_labs)
      
    @eval function Base.unique(ctg::CategoricalArray) 
        l = levels(ctg)
        new_ctg = CategoricalArray(l)
        levels!(new_ctg, l)
    end
  
    plot_months = repeat(month_names, outer = size(m)[2])
    plot_guilds = repeat(guild_labs, inner = size(m)[1])

    plot_mons = categorical(plot_months; levels = month_names)
    plot_glds = categorical(plot_guilds; levels = guild_labs)
  
    return plot_mons, plot_glds
  
end ;

function plot_bar_months(data_matrices, plot_mons, plot_guilds, yl, year)

    lfs = 9
    xtfs = 12
    tfs=16

    p1 = groupedbar(plot_guilds, data_matrices[1], group=plot_mons, 
                    xrotation=45, grid=false, bar_width=0.7, lw=0, 
                    size=(1000,300), bottom_margin=-12mm, left_margin=5mm, top_margin=5mm, xtickfontcolor=:white,
                    legend_position = :bottomright, ylabel="", titlefontsize=tfs, titlelocation=:right,
                    foreground_color_legend = nothing, legendfontsize=lfs, 
                    background_color_legend = nothing, title="5m    ", ylimits=(-40,15));

    p2 = groupedbar(plot_guilds, data_matrices[2], group=plot_mons, 
                    xrotation=45, grid=false, bar_width=0.7, lw=0, left_margin=5mm,
                    size=(1000,300), bottom_margin=-12mm, xtickfontcolor=:white,
                    ylabel="", titlefontsize=tfs, titlelocation=:right,
                    legend=false,  title="DCM    ", ylimits=(-20,12));

    p3 = groupedbar(plot_guilds, data_matrices[3], group=plot_mons, 
                    xrotation=45, grid=false, bar_width=0.7, lw=0, left_margin=5mm,
                    size=(1000,300), bottom_margin=-12mm, xtickfontcolor=:white,
                    ylabel=" % Change in Relative Abundance", ylabelfontsize=16, titlefontsize=tfs, titlelocation=:right,
                    legend=false,  title="150m    ", ylimits=(-6,6));

    p4 = groupedbar(plot_guilds, data_matrices[4], group=plot_mons, 
                    xrotation=45, grid=false, bar_width=0.7, lw=0, left_margin=5mm,
                    size=(1000,300), bottom_margin=-12mm, xtickfontcolor=:white,
                    ylabel="", titlefontsize=tfs, titlelocation=:right,
                    legend=false,  title="500m    ", ylimits=(-6,7));

    p5 = groupedbar(plot_guilds, data_matrices[5], group=plot_mons, 
                    xrotation=45, grid=false, bar_width=0.7, lw=0, 
                    size=(1000,300), bottom_margin=10mm, left_margin=5mm,
                    ylabel="", xlabel="Taxonomic Guild", xlabelfontsize=16, titlelocation=:right,
                    legend=false, titlefontsize=tfs,  title="890m    ", ylimits=(-15,10));

    fig = plot(p1, p2, p3, p4, p5, size=(1000, 800), layout=(5,1), plot_title="Change in RA (2014-2015)", 
    plot_titlefontsize=22, top_margin=-5mm)

    # savefig(fig, "spring_$(year)_guilds.png")
    savefig(fig, "spot_%change_fall_14-15.png")

    return fig

end

#--------------------------------------------------------------------------------------
function get_annual_means(all_in_year, all_clade_colnames)

    clades = all_in_year[!, all_clade_colnames];

    return [mean(skipmissing(c)) for c in eachcol(clades)];

end

function calculate_anomalies(all_in_year, annual_means, all_clade_colnames)

    return all_in_year[!, all_clade_colnames] .- annual_means'

end

function get_anomalies(all_in_year, all_clade_colnames)

    annual_means = get_annual_means(all_in_year, all_clade_colnames);
    anomalies = calculate_anomalies(all_in_year, annual_means, all_clade_colnames);
    anomalies_matrix = create_data_matrices_months(anomalies);

    return anomalies_matrix

end

function plot_monthly_RA_for_year(colnames, matrices, deps, year, copio, slow_copio, oligo)

    dirname = "results/plots/SPOT/monthly_RA_$(year)"
    check_dir_exists(dirname)

    fnames = ["Actinomarinales", "Enterobacterales", "Dadabacteriales","Flavobacteriales","Planctomycetota",
    "Rhodobacterales","SAR11","SAR202","SAR324","SAR406","SAR86","Thioglobaceae","Thermoplasmata","Verrucomicrobiota",
    "Pseudomonadales","Gemmatimonadota","UBA10353","HOC36","Microtrichales","Puniceispirillales", "Rhodospirillales",
    "Ammonia_Ox","Nitrite_Ox"];
    cols=["seagreen3", "seagreen", "darkcyan", "darkblue", "grey20"]
    gnames = ["copio", "slow_copio", "oligo"]
    groups = [copio, slow_copio, oligo]

    a1=0.9
    ylfont=14

    for g in eachindex(groups)

        fsize = size(groups[g])[1]
        n_deps = size(deps)[1]

        f = Vector{Any}(undef, fsize)
        for (i,v) in enumerate(groups[g])
            
            yl = get_bar_limits(v, matrices) 
            dep_plts = Vector{Any}(undef, n_deps)

            for j in eachindex(deps)
                lmar=4mm
                if i == 1
                    dep_plts[j] = bar(matrices[j][:, v], color=cols[j], xlabel="", ylabel=deps[j], left_margin=lmar, xformatter=:none, grid=false, 
                        legend=false, ylimits=yl, lw=0, alpha=a1, ylabelfontsize=ylfont, top_margin=-5mm);
                else
                    dep_plts[j] = bar(matrices[j][:, v], color=cols[j], xlabel="", ylabel="", xformatter=:none, grid=false, 
                        legend=false, ylimits=yl, lw=0, alpha=a1, ylabelfontsize=ylfont);
                end
            end

            tit = fnames[v]
            p = plot(dep_plts..., layout=(n_deps,1), size=(400,500), plot_title="$tit  ")
            f[i] = plot(p)
        
            # check_dir_exists("results/plots/SPOT/monthly_RA_$(year)")
            savefig(p, "results/plots/SPOT/monthly_RA_$(year)/$(fnames[v]).png")
        end

        flength=fsize*250
        fig=plot(f..., layout=(1, fsize), size=(flength, 500))
        savefig(fig, "$dirname/$(gnames[g])_$(year).png")
    end
end


function plot_mean_RA(colnames, matrices, deps, copio, slow_copio, oligo, sems)

    dirname = "results/plots/SPOT/mean_monthly_RA"
    check_dir_exists(dirname)

    fnames = ["Actinomarinales", "Enterobacterales", "Dadabacteriales","Flavobacteriales","Planctomycetota",
    "Rhodobacterales","SAR11","SAR202","SAR324","SAR406","SAR86","Thioglobaceae","Thermoplasmata","Verrucomicrobiota",
    "Pseudomonadales","Gemmatimonadota","UBA10353","HOC36","Microtrichales","Puniceispirillales", "Rhodospirillales",
    "Ammonia_Ox","Nitrite_Ox"];
    cols=["seagreen3", "seagreen", "darkcyan", "darkblue", "grey20"]
    gnames = ["copio", "slow_copio", "oligo"]
    groups = [copio, slow_copio, oligo]

    a1=0.9
    ylfont=14

    for g in eachindex(groups)

        fsize = size(groups[g])[1]
        n_deps = size(deps)[1]

        f = Vector{Any}(undef, fsize)
        for (i,v) in enumerate(groups[g])
            
            yl = get_bar_limits(v, matrices) 
            dep_plts = Vector{Any}(undef, n_deps)

            for j in eachindex(deps)
                lmar=4mm
                if i == 1
                    dep_plts[j] = bar(matrices[j][:, v], color=cols[j], xlabel="", ylabel=deps[j], left_margin=lmar, xformatter=:none, grid=false, 
                        legend=false, ylimits=yl, lw=0, alpha=a1, ylabelfontsize=ylfont, top_margin=-5mm);
                else
                    dep_plts[j] = bar(matrices[j][:, v], color=cols[j], xlabel="", ylabel="", xformatter=:none, grid=false, 
                        legend=false, ylimits=yl, lw=0, alpha=a1, ylabelfontsize=ylfont);
                end
            end

            tit = fnames[v]
            p = plot(dep_plts..., layout=(n_deps,1), size=(400,500), plot_title="$tit  ")
            f[i] = plot(p)
        
            # check_dir_exists("results/plots/SPOT/monthly_RA_$(year)")
            savefig(p, "results/plots/SPOT/mean_monthly_RA/$(fnames[v]).png")
        end

        flength=fsize*250
        fig=plot(f..., layout=(1, fsize), size=(flength, 500))
        savefig(fig, "$dirname/$(gnames[g]).png")
    end
end

function get_bar_limits(i, matrices)

    maxes = Any[];
    for m in matrices
        push!(maxes, unique!([maximum(skipmissing(m[:,i]), init=0) for row in eachrow(m[:, i])])[1])
    end

    maxval = findmax(maxes)[1]
    yl = (0, maxval)
    
    return yl

end

function plot_RA_for_year(months_data, month_nums, all_clade_colnames, copio, slow_copio, oligo, year)

    deps = ["5m", "DCM", "150m", "500m", "890m"];

    dfs = Any[];
    for d in deps
        push!(dfs, get_year_depth_df(months_data, year, d, month_nums))
    end

    matrices = Any[];
    for df in dfs
        push!(matrices, create_data_matrices_months(df, all_clade_colnames))
    end

    plot_monthly_RA_for_year(all_clade_colnames, matrices, deps, year, copio, slow_copio, oligo)

end


function plot_mean_monthly_RA(means_data, month_nums, all_clade_colnames, copio, slow_copio, oligo, sem_data)

    deps = [5, 0, 150, 500, 890];

    means_by_depth = Any[];
    sems_by_depth = Any[];
    for d in deps
        push!(means_by_depth, get_year_depth_df(means_data, year, d, month_nums, true))
        push!(sems_by_depth, get_year_depth_df(sem_data, year, d, month_nums, true))
    end

    matrices = Any[];
    for df in means_by_depth
        push!(matrices, create_data_matrices_months(df, all_clade_colnames))
    end

    plot_mean_RA(all_clade_colnames, matrices, deps, copio, slow_copio, oligo, sems_by_depth)

end



# -----------------------------------------------------------------------------------------------------

year=2015;

month_nums = [1,2,3,4,5,6,7,8,9,10,11,12];
all_clade_colnames = [
"Actinomarinales", "Enterobacterales", "Dadabacteriales","Flavobacteriales","Planctomycetota",
"Rhodobacterales","SAR11_clade","SAR202_clade","SAR324_clade","Marinimicrobia_SAR406_clade",
"SAR86_clade","Thioglobaceae","Thermoplasmata","Verrucomicrobiota","Pseudomonadales",
"Gemmatimonadota","UBA10353_marine_group","HOC36","Microtrichales","Puniceispirillales",
"Rhodospirillales","Ammonia_oxidizers","Nitrite_oxidizers"
];

copio = [2,6,15,5,14,4];
slow_copio = [16,9,19,20,12,10,17,8,13];
oligo = [18,21,11,1,3,7];

monthly_FL_RA = load_monthly_csvs("data/spot_data/monthly_RA/monthly_FL");
plot_RA_for_year(monthly_FL_RA, month_nums, all_clade_colnames, copio, slow_copio, oligo, year);

# MEANS EXCLUDE ANOMALOUS YEAR 2014
# monthly_means_FL_RA, sems = load_monthly_csvs("data/spot_data/monthly_means_RA/", true);

# plot_mean_monthly_RA(monthly_means_FL_RA, month_nums, all_clade_colnames, copio, slow_copio, oligo, sems);









# spring_2014_5 = get_year_depth_df(months_data, 2014, "5m", [3,4,5]);
# spring_2014_dcm = get_year_depth_df(months_data, 2014, "DCM", [3,4,5]);
# spring_2014_150 = get_year_depth_df(months_data, 2014, "150m", [3,4,5]);
# spring_2014_500 = get_year_depth_df(months_data, 2014, "500m", [3,4,5]);
# spring_2014_890 = get_year_depth_df(months_data, 2014, "890m", [3,4,5]);

# spring_2015_5 = get_year_depth_df(months_data, 2015, "5m", [3,4,5]);
# spring_2015_dcm = get_year_depth_df(months_data, 2015, "DCM", [3,4,5]);
# spring_2015_150 = get_year_depth_df(months_data, 2015, "150m", [3,4,5]);
# spring_2015_500 = get_year_depth_df(months_data, 2015, "500m", [3,4,5]);
# spring_2015_890 = get_year_depth_df(months_data, 2015, "890m", [3,4,5]);

# m1 = create_data_matrices_months(spring_2014_5);
# m2 = create_data_matrices_months(spring_2014_dcm);
# m3 = create_data_matrices_months(spring_2014_150);
# m4 = create_data_matrices_months(spring_2014_500);
# m5 = create_data_matrices_months(spring_2014_890);

# m12 = create_data_matrices_months(spring_2015_5);
# m22 = create_data_matrices_months(spring_2015_dcm);
# m32 = create_data_matrices_months(spring_2015_150);
# m42 = create_data_matrices_months(spring_2015_500);
# m52 = create_data_matrices_months(spring_2015_890);

# diff1 = (m12 .- m1) .* 100;
# diff2 = (m22 .- m2) .* 100;
# diff3 = (m32 .- m3) .* 100;
# diff4 = (m42 .- m4) .* 100;
# diff5 = (m52 .- m5) .* 100;

#------------------------------------------------------


# _2014_5 = get_year_depth_df(months_data, 2014, "5m", [9,10,11]);
# _2014_dcm = get_year_depth_df(months_data, 2014, "DCM", [9,10,11]);
# _2014_150 = get_year_depth_df(months_data, 2014, "150m", [9,10,11]);
# _2014_500 = get_year_depth_df(months_data, 2014, "500m", [9,10,11]);
# _2014_890 = get_year_depth_df(months_data, 2014, "890m", [9,10,11]);

# _2015_5 = get_year_depth_df(months_data, 2015, "5m", [9,10,11]);
# _2015_dcm = get_year_depth_df(months_data, 2015, "DCM", [9,10,11]);
# _2015_150 = get_year_depth_df(months_data, 2015, "150m", [9,10,11]);
# _2015_500 = get_year_depth_df(months_data, 2015, "500m", [9,10,11]);
# _2015_890 = get_year_depth_df(months_data, 2015, "890m", [9,10,11]);

# m1 = create_data_matrices_months(_2014_5);
# m2 = create_data_matrices_months(_2014_dcm);
# m3 = create_data_matrices_months(_2014_150);
# m4 = create_data_matrices_months(_2014_500);
# m5 = create_data_matrices_months(_2014_890);

# m12 = create_data_matrices_months(_2015_5);
# m22 = create_data_matrices_months(_2015_dcm);
# m32 = create_data_matrices_months(_2015_150);
# m42 = create_data_matrices_months(_2015_500);
# m52 = create_data_matrices_months(_2015_890);

# diff1 = (m12 .- m1) .* 100;
# diff2 = (m22 .- m2) .* 100;
# diff3 = (m32 .- m3) .* 100;
# diff4 = (m42 .- m4) .* 100;
# diff5 = (m52 .- m5) .* 100;




guild_labs = ["Actinomar", "Enterobac", "Dadabac","Flavobac","Planctomy",
"Rhodobac","SAR11","SAR202","SAR324","SAR406","SAR86",
"Thioglobac","Thermoplas","Verrucomic","Pseudomon","Gemmatim","UBA10353",
"HOC36","Microtric","Puniceispir","Rhodospir","Ammonia_ox","Nitrite_ox"];

# function str_to_datetime(guilds_df)

#     guilds_df.Date = Dates.format.(DateTime.(guilds_df.Date, dateformat"m/d/y") .+ Dates.Year(2000), dateformat"yyyy-mm-dd")
#     guilds_df.Date = Date.(guilds_df.Date)

#     return guilds_df

# end


# function filter_by_depth(guilds_df)

#     guilds_5m = filter(row -> row.Depth == "5m", guilds)
#     guilds_150m = filter(row -> row.Depth == "150m", guilds)

#     return guilds_5m, guilds_150m

# end


# function filter_by_guild()
#     mar2014_sar11 = filter(row -> row.Guild == "SAR11_clade", mar2014)
# end

# path = "data/spot_data/monthly/monthly_FL"
# months_data = load_monthly_csvs(path);

# use Depth.m to split into depths. 
# group rows by month and get monthly averages. Get annual averages too. 
# group by the 4 key guilds emily used 
# plot depths over relative abundance  for annual and seasonal means for each key guild
