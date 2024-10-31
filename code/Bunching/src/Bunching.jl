module Bunching

using DataFrames, CategoricalArrays, Plots, Statistics, CSV, DataFramesMeta, Dates
using LaTeXStrings, Latexify, Random, Printf, Bootstrap, JLD2, AlgebraOfGraphics, RCall

#Script 5: create_sample_map()
function create_sample_map(use_rcall::Bool)    
    #Load all samples and combine
    Fbars = [0.5, 0.6, 0.9, 1.25, 2.0]
    
    df_all = DataFrame()
    for F_bar in Fbars
        fbar_str = replace(string(F_bar), "." => "p")
        println(fbar_str)
        df_tmp = load(joinpath("processed_data", "df_" * fbar_str * ".jld2"),
                               "df_" * fbar_str)
        df_tmp.F_bar .= F_bar
        df_all = vcat(df_all, df_tmp)
    end
    
    dropmissing!(df_all, [:XCoord, :YCoord])
    
    CSV.write(joinpath("processed_data","df_map.csv"), df_all)
    
    ########################################################################################
    ####  Optional: Use RCall to call R script #############################################
    ########################################################################################
    if use_rcall
        R_str = joinpath("code", "scripts", "create_sample_map.R")
        @rput R_str
        R"""
        source(R_str)
        """
    end
    ########################################################################################
end

#Script 4: generate_tables_figures()
function generate_tables_figures(seed::Int)
    ########################################################################################
    #### Set Seed ##########################################################################
    ########################################################################################
    Random.seed!(seed)
    #Base.show(io::IO, f::Float64) = @printf(io, "%.3f", f)
    ########################################################################################
    
    ########################################################################################
    #### Set up directory paths ############################################################
    ########################################################################################
    figures_dir = "figures"
    tables_dir = "tables"
    dpi = 900
    ########################################################################################
    
    ########################################################################################
    #### Produce figures for F̄ = 0.5 #######################################################
    ########################################################################################
    fp_0p5 = joinpath("processed_data", "df_0p5.jld2")
    df_0p5 = load(fp_0p5, "df_0p5")
    δ = 0.04
    F_bar = 0.5
    F_min = 0.1
    F_max = 1.15
    fig_0p5_4delta = Plots.histogram(df_0p5.FAR,  bins = F_min:(δ / 4.0):F_max,
        label = missing, dpi = dpi)
    vline!([F_bar], legend = :none, linewidth = 1.5)
    Plots.ylims!(fig_0p5_4delta, (0, 100))
    xlabel!("F")
    ylabel!("h(F)")
    title!("Figure 7: FAR Distribution for \$\\overline{F}\$ = 0.5")
    fn_svg = "figure_7_FAR_distribution_0p5_4xbins.svg"
    fn_pdf = "figure_7_FAR_distribution_0p5_4xbins.pdf"
    fn_png = "figure_7_FAR_distribution_0p5_4xbins.png"
    fp_svg = joinpath(figures_dir, fn_svg)
    fp_pdf = joinpath(figures_dir, fn_pdf)
    fp_png = joinpath(figures_dir, fn_png)
    savefig(fp_svg)
    savefig(fp_pdf)
    savefig(fp_png)
    ########################################################################################
    
    ########################################################################################
    #### Produce figures for F̄ = 0.6 #######################################################
    ########################################################################################
    fp_0p6 = joinpath("processed_data", "df_0p6.jld2")
    df_0p6 = load(fp_0p6, "df_0p6")
    δ = 0.04
    F_bar = 0.6
    F_min = 0.12
    F_max = 1.5
    fig_0p6_4delta = Plots.histogram(df_0p6.FAR,  bins = F_min:(δ / 4.0):F_max,
        label = missing, dpi = dpi)
    vline!([F_bar], legend = :none, linewidth = 1.5)
    Plots.ylims!(fig_0p6_4delta, (0, 800))
    xlabel!("F")
    ylabel!("h(F)")
    title!("Figure 6: FAR Distribution for \$\\overline{F}\$ = 0.6")
    fn_svg = "figure_6_FAR_distribution_0p6_4xbins.svg"
    fn_pdf = "figure_6_FAR_distribution_0p6_4xbins.pdf"
    fn_png = "figure_6_FAR_distribution_0p6_4xbins.png"
    fp_svg = joinpath(figures_dir, fn_svg)
    fp_pdf = joinpath(figures_dir, fn_pdf)
    fp_png = joinpath(figures_dir, fn_png)
    savefig(fp_svg)
    savefig(fp_pdf)
    savefig(fp_png)
    ########################################################################################
    
    ########################################################################################
    #### Produce figures for F̄ = 0.9 #######################################################
    ########################################################################################
    fp_0p9 = joinpath("processed_data", "df_0p9.jld2")
    df_0p9 = load(fp_0p9, "df_0p9")
    # Plot figures
    δ = 0.07
    F_bar = 0.9
    F_min = F_bar - 10 * δ
    F_max = 2.0
    fig_0p9_4delta = Plots.histogram(df_0p9.FAR,  bins = F_min:(δ / 4.0):F_max,
        label = missing, dpi = dpi)
    vline!([F_bar], legend = :none, linewidth = 1.5)
    Plots.ylims!(fig_0p9_4delta, (0, 215))
    xlabel!("F")
    ylabel!("h(F)")
    title!("Figure 5: FAR Distribution for \$\\overline{F}\$ = 0.9")
    fn_svg = "figure_5_FAR_distribution_0p9_4xbins.svg"
    fn_pdf = "figure_5_FAR_distribution_0p9_4xbins.pdf"
    fn_png = "figure_5_FAR_distribution_0p9_4xbins.png"
    fp_svg = joinpath(figures_dir, fn_svg)
    fp_pdf = joinpath(figures_dir, fn_pdf)
    fp_png = joinpath(figures_dir, fn_png)
    savefig(fp_svg)
    savefig(fp_pdf)
    savefig(fp_png)
    ########################################################################################
    
    ########################################################################################
    #### Produce figures for F̄ = 1.25 ######################################################
    ########################################################################################
    fp_1p25 = joinpath("processed_data", "df_1p25.jld2")
    df_1p25 = load(fp_1p25, "df_1p25")
    # Plot figures
    δ = 0.1
    F_bar = 1.25
    F_min = 0.35
    F_max = 2.75
    fig_1p25_4delta = Plots.histogram(df_1p25.FAR,  bins = F_min:(δ / 4.0):F_max,
        label = missing, dpi = dpi)
    vline!([F_bar], legend = :none, linewidth = 1.5)
    Plots.ylims!(fig_1p25_4delta, (0, 250))
    xlabel!("F")
    ylabel!("h(F)")
    title!("Figure 4: FAR Distribution for \$\\overline{F}\$ = 1.25")
    fn_svg = "figure_4_FAR_distribution_1p25_4xbins.svg"
    fn_pdf = "figure_4_FAR_distribution_1p25_4xbins.pdf"
    fn_png = "figure_4_FAR_distribution_1p25_4xbins.png"
    fp_svg = joinpath(figures_dir, fn_svg)
    fp_pdf = joinpath(figures_dir, fn_pdf)
    fp_png = joinpath(figures_dir, fn_png)
    savefig(fp_svg)
    savefig(fp_pdf)
    savefig(fp_png)
    ########################################################################################
    
    ########################################################################################
    #### Produce figures for F̄ = 2.0 #######################################################
    ########################################################################################
    fp_2p0 = joinpath("processed_data", "df_2p0.jld2")
    df_2p0 = load(fp_2p0, "df_2p0")
    # Plot figures
    δ = 0.15
    F_bar = 2.0
    F_min = 0.35
    F_max = 4.0
    fig_2p0_4delta = Plots.histogram(df_2p0.FAR,  bins = F_min:(δ / 4.0):F_max,
        label = missing, dpi = dpi)
    vline!([F_bar], legend = :none, linewidth = 1.5)
    Plots.ylims!(fig_2p0_4delta, (0, 150))
    xlabel!("F")
    ylabel!("h(F)")
    title!("Figure 3: FAR Distribution for \$\\overline{F}\$ = 2.0")
    fn_svg = "figure_3_FAR_distribution_2p0_4xbins.svg"
    fn_pdf = "figure_3_FAR_distribution_2p0_4xbins.pdf"
    fn_png = "figure_3_FAR_distribution_2p0_4xbins.png"
    fp_svg = joinpath(figures_dir, fn_svg)
    fp_pdf = joinpath(figures_dir, fn_pdf)
    fp_png = joinpath(figures_dir, fn_png)
    savefig(fp_svg)
    savefig(fp_pdf)
    savefig(fp_png)
    ########################################################################################
    
    ########################################################################################
    #### Produce Bunching Illustration Figure: F̄ = 0.6 #####################################
    ########################################################################################
    F_bar = 0.6
    δ = 0.04
    F_min = 0.12
    F_max = 1.5
    FAR_data = df_0p6.FAR[df_0p6.FAR .<= F_max .&& df_0p6.FAR .>= F_min]
    bin_brks = range(F_min, F_max, step = δ)
    p_0p6 = Plots.histogram(FAR_data, alpha = 0.2, bins = bin_brks,
        label = missing, color = :blue, normalize = :pdf, dpi = dpi)
    Plots.ylims!(p_0p6, (0, 4))
    vline!([F_bar - 2*δ, F_bar - δ, F_bar, F_bar + δ, F_bar + 2*δ], label = missing,
        linestyle = :dash, color = :black)
    Plots.xticks!(([F_bar - 2*δ, F_bar, F_bar + 2*δ],
                   [L"\overline{F} - 2 \delta",L"\overline{F}",L"\overline{F} + 2 \delta"]))
    h_min = sum(FAR_data .< (F_bar - δ) .&& FAR_data .>= (F_bar - 2*δ)) /
        size(FAR_data,1) / δ
    h_plus = sum(FAR_data .< (F_bar + 2*δ) .&& FAR_data .>= (F_bar + δ)) /
        size(FAR_data,1)/ δ
    b_1 = sum(FAR_data .< (F_bar) .&& FAR_data .>= (F_bar - δ)) / size(FAR_data,1) / δ
    b_2 = sum(FAR_data .< (F_bar + δ) .&& FAR_data .>= (F_bar)) / size(FAR_data,1) / δ
    #Plots.ylabel!(L"h(F)")
    Plots.yticks!(([h_min, h_plus],[L"h_{-}", ""]))
    Plots.xlims!(p_0p6, (0, F_max))
    annotate!(-0.023,h_min-.13,text(L"h_{+}", :black, :right, 8))
    x_values = [-0, F_bar - δ]
    y_values = [h_min, h_min]
    Plots.plot!(p_0p6, x_values, y_values, color=:green, linestyle = :dash, label = missing)
    x_fill_low = range(F_bar - δ, F_bar, 1000)
    lb_low = [h_min for x in x_fill_low]
    ub_low = [b_1 for x in x_fill_low]
    Plots.plot!(x_fill_low, lb_low, fillrange = ub_low, label = missing, color = :red,
        alpha = 0.5, fillstyle = :/)
    x_fill_high = range(F_bar, F_bar + δ, 1000)
    lb_high = [h_plus for x in x_fill_low]
    ub_high = [b_2 for x in x_fill_low]
    Plots.plot!(x_fill_high, lb_high, fillrange = ub_high, label = L"B", color = :red,
                alpha = 0.5, fillstyle = :/, fillcolor = :red)
    x_values = [-0, F_bar + 2*δ]
    y_values = [h_plus, h_plus]
    Plots.plot!(p_0p6, x_values, y_values, color=:green, linestyle = :dash, label = missing)
    annotate!(1.45, -0.13, L"F")
    Plots.title!("Figure 2: Bunching Area")
    fn_svg = "figure_2_bunching_areas_0p6.svg"
    fn_pdf = "figure_2_bunching_areas_0p6.pdf"
    fn_png = "figure_2_bunching_areas_0p6.png"
    fp_svg = joinpath(figures_dir, fn_svg)
    fp_pdf = joinpath(figures_dir, fn_pdf)
    fp_png = joinpath(figures_dir, fn_png)
    savefig(fp_svg)
    savefig(fp_pdf)
    savefig(fp_png)
    ########################################################################################
    
    ########################################################################################
    #### Produce Table 1: FAR Groups #######################################################
    ########################################################################################
    F_in = [2.0, 1.25, 0.9, 0.6, 0.5]
    obs_in = [size(df_2p0,1), size(df_1p25,1), size(df_0p9,1), size(df_0p6,1),
              size(df_0p5,1)]
    floors_in = [mean(skipmissing(df_2p0.NumFloors)),
                 mean(skipmissing(df_1p25.NumFloors)),
                 mean(skipmissing(df_0p9.NumFloors)),
                 mean(skipmissing(df_0p6.NumFloors)),
                 mean(skipmissing(df_0p5.NumFloors))]
    floor_far_in = [mean(skipmissing(df_2p0.NumFloors ./ df_2p0.FAR)),
                    mean(skipmissing(df_1p25.NumFloors ./ df_1p25.FAR)),
                    mean(skipmissing(df_0p9.NumFloors ./ df_0p9.FAR)),
                    mean(skipmissing(df_0p6.NumFloors ./ df_0p6.FAR)),
                    mean(skipmissing(df_0p5.NumFloors ./ df_0p5.FAR))]
    df_tbl1 = DataFrame(F_bar = F_in, Observations = obs_in,
                        Floors = round.(floors_in; digits = 3),
                        Floors_Over_FAR = round.(floor_far_in; digits = 3))
    # Write to CSV
    fp_tbl1 = joinpath(tables_dir, "table1.csv")
    CSV.write(fp_tbl1, df_tbl1)
    ########################################################################################
    
    ########################################################################################
    #### Produce Table 2: Estimated θ values ###############################################
    ########################################################################################
    F_in = [2.0, 1.25, 0.9, 0.6, 0.5]
    obs_in = [size(df_2p0,1), size(df_1p25,1), size(df_0p9,1), size(df_0p6,1),
              size(df_0p5,1)]
    δ_in = [0.15, 0.1, 0.15, 0.04, 0.04]
    
    #Estimate θ
    N_boot = 10000
    cil = 0.95
    λ = 0.85
    
    θ_2p0, θ_2p0_std, θ_2p0_bias, θ_2p0_bci, θ_2p0_bs =
        bootstrap_θ(df_2p0.FAR, δ_in[1], F_in[1], N_boot, cil)
    θ_1p25, θ_1p25_std, θ_1p25_bias, θ_1p25_bci, θ_1p25_bs =
        bootstrap_θ(df_1p25.FAR, δ_in[2], F_in[2], N_boot, cil)
    θ_0p9, θ_0p9_std, θ_0p9_bias, θ_0p9_bci, θ_0p9_bs =
        bootstrap_θ(df_0p9.FAR, δ_in[3], F_in[3], N_boot, cil)
    θ_0p6, θ_0p6_std, θ_0p6_bias, θ_0p6_bci, θ_0p6_bs =
        bootstrap_θ(df_0p6.FAR, δ_in[4], F_in[4], N_boot, cil)
    θ_0p5, θ_0p5_std, θ_0p5_bias, θ_0p5_bci, θ_0p5_bs =
        bootstrap_θ(df_0p5.FAR, δ_in[5], F_in[5], N_boot, cil)
    
    θ_in = [θ_2p0, θ_1p25, θ_0p9, θ_0p6, θ_0p5]
    ci_lb_in = [θ_2p0_bci[2], θ_1p25_bci[2], θ_0p9_bci[2], θ_0p6_bci[2], θ_0p5_bci[2]]
    ci_ub_in = [θ_2p0_bci[3], θ_1p25_bci[3], θ_0p9_bci[3], θ_0p6_bci[3], θ_0p5_bci[3]]
    avg_θ_in = [θ_2p0_bias, θ_1p25_bias, θ_0p9_bias, θ_0p6_bias, θ_0p5_bias] .+ θ_in
    mc_ratio_in = θ_in .^ λ
    
    bs_mc_ratio = [θ_2p0_bs .^ λ, θ_1p25_bs .^ λ, θ_0p9_bs .^ λ, θ_0p6_bs .^ λ,
                   θ_0p5_bs .^ λ]
    ci_mc_ratio_lb = [quantile(bs_mc_ratio[1], 0.025), quantile(bs_mc_ratio[2], 0.025),
                      quantile(bs_mc_ratio[3], 0.025), quantile(bs_mc_ratio[4], 0.025),
                      quantile(bs_mc_ratio[5], 0.025)]
    ci_mc_ratio_ub = [quantile(bs_mc_ratio[1], 0.975), quantile(bs_mc_ratio[2], 0.975),
                      quantile(bs_mc_ratio[3], 0.975), quantile(bs_mc_ratio[4], 0.975),
                      quantile(bs_mc_ratio[5], 0.975)]
    
    df_tbl2 = DataFrame(F_bar = F_in, Observations = obs_in, delta = δ_in,
        theta_hat = round.(θ_in; digits = 3), ci_lb = round.(ci_lb_in, digits = 3),
        ci_ub = round.(ci_ub_in; digits = 3), avg_theta_hat = round.(avg_θ_in; digits = 3),
        avg_mc_ratio = round.(mc_ratio_in; digits = 3),
        ci_mc_ratio_lb = round.(ci_mc_ratio_lb; digits = 3),
        ci_mc_ratio_ub = round.(ci_mc_ratio_ub; digits = 3))
    #Write to CSV
    CSV.write(joinpath(tables_dir, "table2.csv"), df_tbl2)
    ########################################################################################
    
    ########################################################################################
    #### Produce Table 3: Sensitivity Analysis #############################################
    ########################################################################################

    #Perform sensitivity analaysis
    F_in = [2.0, 2.0, 2.0, 1.25, 1.25, 1.25, 0.9, 0.9, 0.9, 0.6, 0.6, 0.6, 0.5, 0.5, 0.5]
    δ_in = [0.15, 0.1, 0.2, 0.1, 0.05, 0.15, 0.15, 0.125, 0.175, 0.04, 0.03, 0.05, 0.04, 0.03,
            0.05]
    
    θ_2p0_l, ~, ~, ~ = bootstrap_θ(df_2p0.FAR, δ_in[2], F_in[2], N_boot, cil)
    θ_2p0_u, ~, ~, ~ = bootstrap_θ(df_2p0.FAR, δ_in[3], F_in[3], N_boot, cil)
    θ_1p25_l, ~, ~, ~ = bootstrap_θ(df_1p25.FAR, δ_in[5], F_in[5], N_boot, cil)
    θ_1p25_u, ~, ~, ~ = bootstrap_θ(df_1p25.FAR, δ_in[6], F_in[6], N_boot, cil)
    θ_0p9_l, ~, ~, ~ = bootstrap_θ(df_0p9.FAR, δ_in[8], F_in[8], N_boot, cil)
    θ_0p9_u, ~, ~, ~ = bootstrap_θ(df_0p9.FAR, δ_in[9], F_in[9], N_boot, cil)
    θ_0p6_l, ~, ~, ~ = bootstrap_θ(df_0p6.FAR, δ_in[11], F_in[11], N_boot, cil)
    θ_0p6_u, ~, ~, ~ = bootstrap_θ(df_0p6.FAR, δ_in[12], F_in[12], N_boot, cil)
    θ_0p5_l, ~, ~, ~ = bootstrap_θ(df_0p5.FAR, δ_in[14], F_in[14], N_boot, cil)
    θ_0p5_u, ~, ~, ~ = bootstrap_θ(df_0p5.FAR, δ_in[15], F_in[15], N_boot, cil)
    
    θ_in = [θ_2p0, θ_2p0_l, θ_2p0_u, θ_1p25, θ_1p25_l, θ_1p25_u, θ_0p9, θ_0p9_l, θ_0p9_u,
            θ_0p6,
            θ_0p6_l, θ_0p6_u, θ_0p5, θ_0p5_l, θ_0p5_u]
    df_tbl3 = DataFrame(F_bar = F_in, delta = δ_in, theta = round.(θ_in; digits = 3))
    CSV.write(joinpath(tables_dir, "table3.csv"), df_tbl3)
    ########################################################################################
    
    ########################################################################################
    #### Table 4: Floor Space Counterfactual ###############################################
    ########################################################################################
    
    #Set-up vectors for table columns
    df_in = [df_2p0, df_1p25, df_0p9, df_0p6, df_0p5]
    F_in = [2.0, 1.25, 0.9, 0.6, 0.5]
    δ_in = [0.15, 0.1, 0.15, 0.04, 0.04]
    θ_in = [θ_2p0, θ_1p25, θ_0p9, θ_0p6, θ_0p5]
    
    #Compute dummies for bunchers, below bunchers, and above bunchers
    for i ∈ eachindex(df_in)
        #Create buncher dummy: F ∈ [F̄ - δ, F̄ + δ]
        df_in[i].buncher .= false
        buncher_cond = (round.(df_in[i].FAR, digits = 2) .<=
            round.(F_in[i] .+ δ_in[i],digits = 2)) .& (round.(df_in[i].FAR, digits = 2) .>=
            (round.(F_in[i] .- δ_in[i], digits = 2)))

        df_in[i][buncher_cond, :buncher] .= true
    
        #Create below and above buncher dummies: F < F̄ - δ
        df_in[i].below_buncher .= false
        below_buncher_cond = (round.(df_in[i].FAR, digits = 2) .<
            round(F_in[i] .- δ_in[i], digits = 2))
        df_in[i][below_buncher_cond, :below_buncher] .= true
    
        #Create above buncher dummy: F > F̄ + δ
        df_in[i].above_buncher .= false
        above_buncher_cond = (round.(df_in[i].FAR, digits = 2) .>
            round(F_in[i] .+ δ_in[i], digits = 2))
        df_in[i][above_buncher_cond, :above_buncher] .= true
    
        #Compute spacenow
        df_in[i].spacenow .= 0.0
        df_in[i][df_in[i].below_buncher, :spacenow] .=
            df_in[i][df_in[i].below_buncher, :].LotArea .*
            df_in[i][df_in[i].below_buncher, :].FAR
        df_in[i][df_in[i].above_buncher, :spacenow] .=
            df_in[i][df_in[i].above_buncher, :].LotArea .*
            df_in[i][df_in[i].above_buncher, :].FAR
        df_in[i][df_in[i].buncher, :spacenow] .=
            df_in[i][df_in[i].buncher, :].LotArea .* F_in[i]
        #if any spacenow_in == 0.0 throw error
        if any(df_in[i].spacenow .== 0.0)
            println("Error: spacenow_in is zero for F̄ = ", F_in[i])
        end
    
        #Compute spacenew
        df_in[i].spacenew .= 0.0
        df_in[i][df_in[i].below_buncher, :spacenew] .= 
            df_in[i][df_in[i].below_buncher, :spacenow]
        df_in[i][df_in[i].above_buncher, :spacenew] .=
            df_in[i][df_in[i].above_buncher, :].LotArea .*
            df_in[i][df_in[i].above_buncher, :].FAR .* θ_in[i]
        df_in[i][df_in[i].buncher, :spacenew] .=
            0.5 .* df_in[i][df_in[i].buncher, :].LotArea .* (1.0 .+ θ_in[i]) .* F_in[i]
        #if any spacenew_in == 0.0 throw error
        if any(df_in[i].spacenew .== 0.0)
            println("Error: spacenow_in is zero for F̄ = ", F_in[i])
        end
    
    end
    
    spacenow_in = [sum(skipmissing(df_in[i].spacenow)) for i ∈ eachindex(df_in)]
    spacenew_in = [sum(skipmissing(df_in[i].spacenew)) for i ∈ eachindex(df_in)]
    percent_gain_in = (spacenew_in .- spacenow_in) ./ spacenow_in
    
    
    
    df_tbl4 = DataFrame(F_bar = F_in, delta = δ_in, Space_Now = round.(spacenow_in; digits = 3),
        Space_New = round.(spacenew_in; digits = 3),
        percent_gain = round.(percent_gain_in; digits = 3))
    
    #Write to CSV
    fp_tbl4 = joinpath(tables_dir, "table4.csv")
    CSV.write(fp_tbl4, df_tbl4)
    ########################################################################################
    
    ########################################################################################
    #### Produce Figure 9: Bootstrap Distribution of θ for F̄ = 2.0 #########################
    ########################################################################################
    #Compute Bootstrap Historgram
    f_bs_2p0 = Plots.histogram(θ_2p0_bs, label = missing, dpi = dpi)
    xlabel!(L"\hat{\theta}")
    vline!([θ_2p0], label = "Point Estimate", linestyle = :dash, linewidth = 2)
    Plots.ylims!(f_bs_2p0, (0, 900))
    Plots.xlims!(f_bs_2p0, (1.045, 1.225))
    
    title!("Figure 9: Bootstrap Distribution of \$\\hat{\\theta}\$ for " *
           "\$\\overline{F}\$ = 2.0")
    fn_svg = "figure_9_bootstrap_hist_2p0.svg"
    fn_pdf = "figure_9_bootstrap_hist_2p0.pdf"
    fn_png = "figure_9_bootstrap_hist_2p0.png"
    savefig(joinpath(figures_dir, fn_svg))
    savefig(joinpath(figures_dir, fn_pdf))
    savefig(joinpath(figures_dir, fn_png))
    ########################################################################################
end

#Script 3: generate_analysis_subsamples()
function generate_analysis_subsamples(df::DataFrame) 
    ########################################################################################
    #### Load Required DataFrames ##########################################################
    ########################################################################################
    #if !(@isdefined df)
    #    fp_df_cleaned = joinpath("processed_data", "df_transformed_cleaned_pluto.jld2")
    #    df = load(fp_df_cleaned, "df")
    #end
    ########################################################################################
    
    ########################################################################################
    #### FBar = 0.5 ########################################################################
    ########################################################################################
    df.YQ = df.Year #Rename Year -> YQ
    #df_fc.YQ = df_fc.Year #Rename Year -> YQ
    
    # Get subsample of data for F̄ = 0.5
    zonestring_0p5 = ["R1-1", "R1-2", "R1-2A","R2", "R2A"] #Zoning districts for F̄ = 0.5
    df_0p5 = generate_subsample(zonestring_0p5, df) #Generate subsample
    
    # Save data
    fp_0p5_jld = joinpath("processed_data", "df_0p5.jld2")
    jldsave(fp_0p5_jld; df_0p5)
    fp_0p5_csv = joinpath("processed_data", "df_0p5.csv")
    CSV.write(fp_0p5_csv, df_0p5)
    ########################################################################################
    
    ########################################################################################
    #### F̄ = 0.6 ###########################################################################
    ########################################################################################
    
    # Get subsample of data for F̄ = 0.6
    zonestring_0p6 = ["R3-1", "R3-2", "R3A","R3X"] #Zoning districts for F̄ = 0.6
    df_0p6 = generate_subsample(zonestring_0p6, df)
    
    # Save Data
    fp_0p6_jld = joinpath("processed_data", "df_0p6.jld2")
    jldsave(fp_0p6_jld; df_0p6)
    fp_0p6_csv = joinpath("processed_data", "df_0p6.csv")
    CSV.write(fp_0p6_csv, df_0p6)
    ########################################################################################
    
    ########################################################################################
    #### F̄ = 0.9 ###########################################################################
    ########################################################################################
    zonestring_0p9 = ["R4", "R4-1", "R4A", "R4B"]
    df_0p9 = generate_subsample(zonestring_0p9, df)
    fp_0p9_csv = joinpath("processed_data", "df_0p9.csv")
    CSV.write(fp_0p9_csv, df_0p9)
    fp_0p9_jld = joinpath("processed_data", "df_0p9.jld2")
    jldsave(fp_0p9_jld; df_0p9)
    ########################################################################################
    
    ########################################################################################
    #### F̄ = 1.25 ##########################################################################
    ########################################################################################
    zonestring_1p25 = ["R5"]
    df_1p25 = generate_subsample(zonestring_1p25, df)
    fp_1p25_csv = joinpath("processed_data", "df_1p25.csv")
    CSV.write(fp_1p25_csv, df_1p25)
    fp_1p25_jld = joinpath("processed_data", "df_1p25.jld2")
    jldsave(fp_1p25_jld; df_1p25)
    ########################################################################################
    
    ########################################################################################
    #### F̄ = 2.0 ###########################################################################
    ########################################################################################
    zonestring_2p0 = ["R5D", "R6B"]
    df_2p0 = generate_subsample(zonestring_2p0, df)
    fp_2p0_csv = joinpath("processed_data", "df_2p0.csv")
    CSV.write(fp_2p0_csv, df_2p0)
    fp_2p0_jld = joinpath("processed_data", "df_2p0.jld2")
    jldsave(fp_2p0_jld; df_2p0)
    ########################################################################################
    
end

#Script 2: transform_cleaned_pluto_data()
function transform_cleaned_pluto_data(df::DataFrame, save_results::Bool)

    table_dir = joinpath("tables")
    
    # Load the data
    #if !(@isdefined df)
    #    fp_pluto_clean = joinpath("processed_data", "df_cleaned_pluto.jld2")
    #    df = jldopen(fp_pluto_clean) do file
    #        file["df_pluto"]
    #    end
    #end
    
    
    #Save size of initial dataset
    N_1 = nrow(df[df.Year .== 2017,:]) #N_1 = 859,223 
    
    #Get intital sample size
    boro_counts = combine(groupby(df[df.Year .== 2017,:], [:BoroCode]), nrow => :Count)
    boro_counts.Percentages = boro_counts.Count ./ N_1
    fp_boro_counts = joinpath(table_dir, "initial_boro_counts.csv")
    CSV.write(fp_boro_counts, boro_counts)
    
    #Define CondoDum
    transform!(df, :BldgClass => ByRow(def_condodum) => :CondoDum)
    
    #Drop if YearBuilt is missing, not equal to "E", and not equal to "                  @                     "
    filter!(row -> !ismissing(row.YearBuilt), df)
    filter!(row -> row.YearBuilt != "E", df)
    filter!(row -> row.YearBuilt != "                  @                     ", df)
    
    #Correctly Parse YearBuilt to integer
    transform!(df, :YearBuilt => ByRow(correct_year_built) => :YearBuilt)
    
    # How many properties increased height / FAR in the dataset?
    # Removed YQ as it is not present in new dataset
    select!(df, :BoroCode, :Block, :Lot, :Address, :CT2010, :CB2010, :NumBldgs,
        :NumFloors, :BldgArea, :BldgDepth, :BldgFront, :BldgClass, :LotArea,
        :LotDepth, :LotFront, :ComArea, :ResArea, :OfficeArea, :RetailArea,
        :GarageArea, :StrgeArea, :FactryArea, :UnitsRes, :UnitsTotal,
        :YearBuilt, :YearAlter1, :YearAlter2, :IrrLotCode, :LandUse,
        :FAR, :MaxAllwFAR, :ResidFAR, :CommFAR, :FacilFAR, :CornerLot,
        :CondoNo, :LotType, :LtdHeight, :ZoneDist1, :ZoneDist2, :ZoneDist3,
        :ZoneDist4, :Overlay1, :Overlay2, :SPDist1, :SPDist2, :SPDist3,
        :SplitZone, :HistDist, :CondoDum, :Year, :XCoord, :YCoord, :PlutoVer)
        
    #Drop Split Zones
    select!(df, Not(:SplitZone))
    filter!(row -> ismissing(row.Overlay2), df)
    select!(df, Not(:Overlay2))
    filter!(row -> ismissing(row.SPDist2), df)
    select!(df, Not(:SPDist2))
    filter!(row -> ismissing(row.SPDist3), df)
    select!(df, Not(:SPDist3))
    
    #Drop any special purpose districts
    filter!(row -> ismissing(row.SPDist1), df)
    select!(df, Not(:SPDist1))
    
    #Drop Historical districts
    filter!(row -> ismissing(row.HistDist), df)
    select!(df, Not(:HistDist))
    
    #Drop Limited Height districts
    filter!(row -> ismissing(row.LtdHeight), df)
    select!(df, Not(:LtdHeight))
        
    # Filter on Num Bldgs<2, NumFloors > 0,
    filter!(row -> !(ismissing(row.NumBldgs)), df)
    filter!(:NumBldgs => <(2), df)
    
    #See what percent of observations have been lost thus far.
    N_2 = nrow(df[df.Year .== 2017, :]) #N_2 = 484,021 
    diff_1 = N_1 - N_2 #diff_1 = 375,202
    per_loss_1 = diff_1 / N_1 #per_loss_1 = 0.437
    
    # Filter BldgClass not in {G*, I*, J*, L*, M*, N*, P*, Q*, T*, U*, W*, Y*, Z*}
    filter!(row -> !(ismissing(row.BldgClass)), df)
    bldg_class_filter = ["G", "I", "J", "L", "M", "N", "P", "Q", "T", "U", "W", "Y", "Z"]
    tmp_fn(x) = filter_bldg_class(x, bldg_class_filter)
    filter!(:BldgClass .=> tmp_fn, df)
    
    #Drop irregular Zones
    dropmissing!(df, :ZoneDist1)
    filter!(row -> row.ZoneDist1 != "BPC", df)
    #df = df[Not(df.ZoneDist1 .== "BPC"),:]
    filter!(row -> row.ZoneDist1 != "NZS", df)
    #df = df[Not(df.ZoneDist1 .== "NZS"),:]
    filter!(row -> row.ZoneDist1 != "PARK", df)
    #df = df[Not(df.ZoneDist1 .== "PARK"),:]
    filter!(row -> row.ZoneDist1 != "PARKNYS", df)
    #df = df[Not(df.ZoneDist1 .== "PARKNYS"),:]
    filter!(row -> row.ZoneDist1 != "PARKNY", df)
    filter!(row -> row.ZoneDist1 != "Drop Lot", df)
    filter!(row -> row.ZoneDist1 != "PARKUS", df)
    
    
    # Your conversion to Vector{String} is fine, but it's not typically necessary unless you're ensuring type stability for some specific reason.
    df.ZoneDist1 = Vector{String}(df.ZoneDist1)
    
    # Use filter! to modify df in place, removing the rows where ZoneDist1 contains "/"
    filter!(row -> is_mixed_zone(row.ZoneDist1), df)
    
    #Drop commercial
    filter!(row -> string(row.ZoneDist1[1]) != "C", df)
    filter!(row -> string(row.ZoneDist1[1]) != "M", df)
    
    N_3 = nrow(df[df.Year .== 2017, :]) #N_3 = 445,302
    
    #Generate attic and IH allowances
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "M1-1D" => "M1-1")) => :ZoneDist1)
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "M1-2D" => "M1-2")) => :ZoneDist1)
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "M1-3D" => "M1-3")) => :ZoneDist1)
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "M1-4D" => "M1-4")) => :ZoneDist1)
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "M1-5A" => "M1-5")) => :ZoneDist1)
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "M1-5B" => "M1-5")) => :ZoneDist1)
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "M1-5M" => "M1-5")) => :ZoneDist1)
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "M1-6M" => "M1-6")) => :ZoneDist1)
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "M1-6D" => "M1-6")) => :ZoneDist1)
    transform!(df, :ZoneDist1 => ByRow(x -> replace(x, "R10H" => "R10")) => :ZoneDist1)
    filter!(:ZoneDist1 => !=("R7-3"), df)
    filter!(:ZoneDist1 => !=("ZR 11-151"), df)
    
    transform!(df, :ZoneDist1 => ByRow(generate_fbar) => :FAR_bar)
    transform!(df, :ZoneDist1 => ByRow(generate_attic_allowance) => :AtticAllowance)
    transform!(df, :ZoneDist1 => ByRow(generate_ih_bonus) => :IHBonus)
    transform!(df, [:ZoneDist1, :Overlay1] => ByRow(generate_overlay_bonus) => :OverlayBonus)
    transform!(df, :ZoneDist1 => ByRow(generate_community_facility_fbar) => :ComFacilFAR)
    
    #If zoned M, and have an overlay, drop row
    filter!([:ZoneDist1, :Overlay1] => drop_manufacturing_overlay, df)
    transform!(df, [:FAR_bar, :AtticAllowance] => ((x, y) -> x + y) => :FAR_bar_attic)
    transform!(df, [:FAR_bar, :IHBonus] => ((x, y) -> x + y) => :FAR_bar_ih)
    transform!(df, [:FAR_bar, :OverlayBonus] => ((x, y) -> x + y) => :FAR_bar_overlay)
    transform!(df, [:FAR_bar, :AtticAllowance, :OverlayBonus] => ((x, y, z) -> x + y + z) => :FAR_bar_attic_overlay)
    transform!(df, [:FAR_bar, :IHBonus, :OverlayBonus] => ((x, y, z) -> x + y + z) => :FAR_bar_ih_overlay)
    transform!(df, [:FAR_bar, :AtticAllowance, :ComFacilFAR] => ((x, y, z) -> x + y + z) => :FAR_bar_attic_com)
    transform!(df, [:FAR_bar, :IHBonus, :ComFacilFAR] => ((x, y, z) -> x + y + z) => :FAR_bar_ih_com)
    transform!(df, [:FAR_bar, :OverlayBonus, :ComFacilFAR] => ((x, y, z) -> x + y + z) => :FAR_bar_overlay_com)
    transform!(df, [:FAR_bar, :IHBonus, :OverlayBonus, :ComFacilFAR] => ((x, y, z, w) -> x + y + z + w) => :FAR_bar_ih_overlay_com)
    
    #If MaxAllwFAR is missing, set it to FAR with majority LandUse
    transform!(df, :MaxAllwFAR => ByRow(correct_max_allw_far) => :MaxAllwFAR)
    
    transform!(df, 
        [:MaxAllwFAR, :ComArea, :ResArea, :ResidFAR, :CommFAR] => 
        ByRow(xform_maxfar) => :MaxAllwFAR)
    
    #Drop missing values of MaxAllwFAR
    dropmissing!(df, :MaxAllwFAR)
    
    # Filter on MaxAllwFAR > 0
    filter!(:MaxAllwFAR => >(0.), df)
    
    # Round Down NumFloors to nearest integer
    df[!, :NumFloors] = floor.(df[!, :NumFloors])
    
    #For each building, set numfloors to maximum number of floors
    transform!(groupby(df, [:BoroCode,:Block,:Lot,:BldgArea]),
        :NumFloors => maximum => :NumFloors)
    
    transform!(df, :YearAlter1 => ByRow(correct_year_built) => :YearAlter1)
    transform!(df, :YearAlter2 => ByRow(correct_year_built) => :YearAlter2)
    
    
    #Create Ratio of FAR / MaxAllwFAR
    dropmissing!(df, :FAR)
    
    #Convert BoroCode to Int64
    transform!(df, :BoroCode=> ByRow(correct_boro_code) => :BoroCode)
    transform!(df, :Block => ByRow(correct_boro_code) => :Block)
    transform!(df, :Lot => ByRow(correct_boro_code) => :Lot)
    
    
    #Get SubDataFrame with Year == 2017 and YearBuilt >= 2000
    df_y17 = @rsubset df begin
        :Year == 2017
        :YearBuilt >= 2000
    end
    
    tmp_fn(BoroCode::Int64, Block::Int64, Lot::Int64, YearBuilt::Int64) = 
        zoning_dist_year_built(BoroCode, Block, Lot, YearBuilt, df)
    
    # Apply the function row-wise
    @rtransform! df_y17 :ZoneDistYearBuilt = tmp_fn(:BoroCode, :Block, :Lot, :YearBuilt)

    #Save processed cleaned data
    df = nothing
    df = df_y17
    
    if save_results
        fp_df_cleaned = joinpath("processed_data", "df_transformed_cleaned_pluto.jld2")
        jldsave(fp_df_cleaned; df)
    end

    return df
end

#Script 1: clean_raw_pluto_data.jl
function clean_raw_pluto_data(save_results::Bool)

    pluto_root = joinpath("raw_data", "PLUTO")
    
    ######    Load in metadata on column names throughout    #######################
    fp_vc = joinpath(pluto_root, "extra_csv_files", "variableChart.csv")
    variable_chart = CSV.read(fp_vc, DataFrame)
    
    varnames_in = names(variable_chart)[3:end]
    varnames_out = vcat(varnames_in, ["Year", "PlutoVer"])
    
    fp_ybf = joinpath(pluto_root, "extra_csv_files", "year_borocode_fp.csv")
    year_bc_fp = CSV.read(fp_ybf, DataFrame)
    
    fp_ar = joinpath(pluto_root, "extra_csv_files", "actionRequired.csv")
    action_required = CSV.read(fp_ar, DataFrame)
    
    #####    Main Loop    #####
    create_df_trig = 0
    year_init = 2002
    df_pluto = DataFrame()
    year_bc_fp = year_bc_fp[year_bc_fp.Year .>= year_init,:]
    for (n, (yr, boro_code, fp, pluto_ver)) in enumerate(eachrow(year_bc_fp))
        println("Year: $yr, BoroCode: $boro_code, PlutoVersion: $pluto_ver")
        #Get file path
        fp = year_bc_fp[
            (year_bc_fp.Year .== yr) .& (year_bc_fp.BoroCode .== boro_code),
            :FilePath][1]
    
        #Load data
        df_tmp = CSV.read(joinpath(pluto_root, fp), DataFrame)
    
        df_tmp.Year .= yr
        df_tmp.PlutoVer .= pluto_ver
    
        #Initialize special_case triggers
        special_2_trig = 0
        special_3_trig = 0
    
        #Loop over columns and apply transformations
        for name_out in varnames_in
            
            println("Variable Name: $name_out")
            name_in =
                variable_chart[variable_chart.PlutoID .== pluto_ver,name_out][1]
    
            #Skip if name_in is "0"
            if name_in == "0"
                df_tmp[!,Symbol(name_out)] .= missing
                continue
            end
    
            #If not, then apply appropiate treansformation
            transform_type =
                action_required[variable_chart.PlutoID .== pluto_ver, name_out][1]
    
            #Skip if tranform_type == "N"
            if transform_type == "N"
                #If name_out != name_in, transform name_in to name_out
                if name_out != name_in
                    transform!(df_tmp, Symbol(name_in) => Symbol(name_out))
                end
                continue
            end
    
            #If not, then apply appropiate treansformation
            if transform_type == "strip"
                transform!(df_tmp, Symbol(name_in) => 
                    ByRow(x -> strip_missing(x) == "" ? missing : strip(x)) => 
                    Symbol(name_out))
    
            elseif transform_type == "strip-Int64"
                transform!(df_tmp, Symbol(name_in) => 
                    ByRow(x -> filter(isdigit,strip_missing(string(x))) == "" ? missing :
                    parse(Int64, filter(isdigit,strip(string(x))))) => 
                    Symbol(name_out))
    
            elseif transform_type == "Special-1"
                transform!(df_tmp, [:NumBldgs, Symbol(name_in)] => 
                    ByRow(special_1) => Symbol(name_out))
    
            elseif transform_type == "zero"
                transform!(df_tmp, Symbol(name_in) => 
                    ByRow(x ->  ismissing(x) ? missing : ifelse(x==0,missing,x)) => 
                    Symbol(name_out))
    
            elseif transform_type == "Special-2"
                if special_2_trig == 0
                    name_1 = variable_chart[variable_chart.PlutoID .== pluto_ver,
                        "ComArea"][1]
    
                    name_2 = variable_chart[variable_chart.PlutoID .== pluto_ver,
                        "ResArea"][1]
    
                    transform!(df_tmp,
                        [:NumBldgs, Symbol(name_1), Symbol(name_2)] => 
                        ByRow(special_2) => [:ComArea, :ResArea])
                    
                    special_2_trig = 1
                end
    
            elseif transform_type == "Special-3"
                if special_3_trig == 0
                    name_1 = variable_chart[variable_chart.PlutoID .== pluto_ver,
                        "ComArea"][1]
                    name_2 = variable_chart[variable_chart.PlutoID .== pluto_ver,
                        "ResArea"][1]
                    name_3 = variable_chart[variable_chart.PlutoID .== pluto_ver,
                        "FactryArea"][1]
                    name_4 = variable_chart[variable_chart.PlutoID .== pluto_ver,
                        "GarageArea"][1]
                    name_5 = variable_chart[variable_chart.PlutoID .== pluto_ver,
                        "RetailArea"][1]
                    name_6 = variable_chart[variable_chart.PlutoID .== pluto_ver,
                        "StrgeArea"][1]
                    name_7 = variable_chart[variable_chart.PlutoID .== pluto_ver,
                        "OfficeArea"][1]
    
                    #Parse each as a Int
                    transform!(df_tmp, Symbol(name_1) => 
                        ByRow(x -> filter(isdigit,strip_missing(x))  == "" ? missing :
                        parse(Int64, filter(isdigit,strip(string(x))))) => 
                        Symbol(name_1))
                    transform!(df_tmp, Symbol(name_2) => 
                        ByRow(x -> filter(isdigit,strip_missing(x))  == "" ? missing :
                        parse(Int64, filter(isdigit,strip(string(x))))) => 
                        Symbol(name_2))
                    transform!(df_tmp, Symbol(name_3) =>
                        ByRow(x -> filter(isdigit,strip_missing(x))  == "" ? missing :
                        parse(Int64, filter(isdigit,strip(string(x))))) => 
                        Symbol(name_3))
                    transform!(df_tmp, Symbol(name_4) =>
                        ByRow(x -> filter(isdigit,strip_missing(x))  == "" ? missing :
                        parse(Int64, filter(isdigit,strip(string(x))))) => 
                        Symbol(name_4))
                    transform!(df_tmp, Symbol(name_5) =>
                        ByRow(x -> filter(isdigit,strip_missing(x))  == "" ? missing :
                        parse(Int64, filter(isdigit,strip(string(x))))) => 
                        Symbol(name_5))
                    transform!(df_tmp, Symbol(name_6) =>
                        ByRow(x -> filter(isdigit,strip_missing(x))  == "" ? missing :
                        parse(Int64, filter(isdigit,strip(string(x))))) => 
                        Symbol(name_6))
                    transform!(df_tmp, Symbol(name_7) =>
                        ByRow(x -> filter(isdigit,strip_missing(x))  == "" ? missing :
                        parse(Int64, filter(isdigit,strip(string(x))))) => 
                        Symbol(name_7))
                    
    
                    transform!(df_tmp, [:NumBldgs, Symbol(name_1), Symbol(name_2),
                        Symbol(name_3), Symbol(name_4), Symbol(name_5),
                        Symbol(name_6), Symbol(name_7)] => 
                        ByRow(special_3) => [:ComArea, :ResArea, :FactryArea,
                            :GarageArea, :RetailArea, :StrgeArea, :OfficeArea])
    
                    special_3_trig = 1
                end
    
            elseif transform_type == "strip-YN"
                transform!(df_tmp, Symbol(name_in) => 
                    ByRow(x -> strip_missing(x) == "" ? missing :
                        ifelse(strip(x) == "Y", true, false)) =>
                    Symbol(name_out))
    
            elseif transform_type == "clean_zone_dist"
                transform!(df_tmp, Symbol(name_in) => 
                    ByRow(clean_zone_dist) => Symbol(name_out))
                
            elseif transform_type == "strip-filterstr-Int64"
                transform!(df_tmp, Symbol(name_in) => 
                    ByRow(x -> filter(isdigit,strip_missing(string(x)))  == "" ? missing :
                    parse(Int64, filter(isdigit,strip(string(x))))) => 
                    Symbol(name_out))
            
            elseif transform_type == "strip-Int64-zero"
                transform!(df_tmp, Symbol(name_in) =>
                    ByRow(x -> 
                        if typeof(x) == String
                            replace(x, "\x10" => "")
                        else
                            x
                        end) => Symbol(name_in))
    
                transform!(df_tmp, Symbol(name_in) => 
                    ByRow(x -> filter(isdigit,strip_missing(string(x)))  == "" ? missing :
                        parse(Int64, filter(isdigit, strip(string(x))))) => 
                    Symbol(name_out))
    
                transform!(df_tmp, Symbol(name_out) => 
                    ByRow(x ->  ismissing(x) ? missing : ifelse(x==0,missing,x)) => 
                    Symbol(name_out))
            else 
                println("Error: $transform_type not recognized")
            end
        end
    
        #Only keep needed columns
        df_tmp = df_tmp[:,map(Symbol,varnames_out)]
    
        #Concatenate dataframes
        if create_df_trig == 0
                df_pluto = df_tmp
                create_df_trig = 1
        else
                df_pluto = vcat(df_pluto, df_tmp)
        end
    end
    
    #Save data
    if save_results
        fp_pluto_clean = joinpath("processed_data", "df_cleaned_pluto.jld2")
        jldsave(fp_pluto_clean; df_pluto)
    end

    return df_pluto
    
end

function clean_zone_dist(zd)
    if ismissing(zd)
        return missing
    end
    szp = strip(zd)
    missing_list = ["", "ZNA", "NZS","DROP LOT"]
    if szp in missing_list
        return missing
    else
        return szp
    end
end

#Special-3
function special_3(nb, ca, ra, fa, ga, reta, stra, offa) 
    if ismissing(nb) | ismissing(ca) | ismissing(ra) | ismissing(fa) |
        ismissing(ga) | ismissing(reta) | ismissing(stra) | ismissing(offa)
        return missing, missing, missing, missing, missing, missing, missing
    elseif nb > 1 && (ca + ra + fa + ga + reta + stra + offa) ≈ 0
        return missing, missing, missing, missing, missing, missing, missing
    else
        return ca, ra, fa, ga, reta, stra, offa
    end
end

function special_2(nb, ca, ra) 

    if ismissing(nb) | ismissing(ca) | ismissing(ra)
        return missing, missing
    elseif nb > 1 && (ca + ra) ≈ 0
        return missing, missing
    else
        return ca, ra
    end
end

function special_1(nb, x) 
    if ismissing(nb) | ismissing(x)
        return missing
    elseif nb > 0. && x == 0.
        return missing
    else
        return x
    end
end

function strip_missing(x)
    if ismissing(x)
        return ""
    else
        return strip(string(x))
    end
end


function count_digits_after_decimal(f::Float64)
    # Convert the float to a string
    s = string(f)
    
    # Find the position of the decimal point
    dot_pos = findfirst('.', s)
    
    # If no decimal point is found, return 0
    if dot_pos === nothing
        return 0
    end
    
    # Return the number of characters after the decimal point
    return length(s) - dot_pos
end


function est_θ(F::Vector{Float64}, δ::Float64, F_bar::Float64)

    #LABELS OBSERVATIONS AS BEING IN minus, plus, OR b (NEAR FBAR) RANGES
    sigdigs = maximum([count_digits_after_decimal(F_bar),
        count_digits_after_decimal(δ)])
    Hminus_dum = (F .<= round(F_bar - δ; digits = sigdigs)) .& (F .> round(F_bar - 2*δ; digits = sigdigs))
    Hplus_dum = (F .>= round(F_bar + δ; digits = sigdigs)) .& (F .< round(F_bar + 2*δ; digits = sigdigs))
    B_dum = (F .> round(F_bar - δ; digits = sigdigs)) .& (F .< round(F_bar + δ; digits = sigdigs))

    #Get sums
    Hminus = sum(Hminus_dum)
    Hplus = sum(Hplus_dum)
    B_pre = sum(B_dum)

    #Get density heights by dividing by d
    Hminus_height = Hminus / δ
    Hplus_height = Hplus / δ

    #GENERATES B AND VARIABLE X APPEARING IN QUADRATIC
    B = B_pre - Hminus - Hplus
    #X = (F_bar*(Hplus_height - Hminus_height) / 2.0) - B
    #Solve for X using Quadratic Formula
    #θ = (-X + (X^2 + F_bar * Hminus_height * Hplus_height)^0.5) / (2.0 * F_bar * Hminus_height)

    a = F_bar * Hminus_height / 2.0
    b = F_bar * Hplus_height / 2.0 - F_bar * Hminus_height / 2.0 - B
    c = -F_bar * Hplus_height / 2.0

    θ = (-b + sqrt(b^2 - 4.0*a*c)) / (2.0 * a)

    #D = Hplus_height * ((1 + θ) / θ) * (θ - 1.0) * F_bar / 2.0

    return θ

end


function est_θ_D(F::Vector{Float64}, δ::Float64, F_bar::Float64)

    #LABELS OBSERVATIONS AS BEING IN minus, plus, OR b (NEAR FBAR) RANGES
    sigdigs = maximum([count_digits_after_decimal(F_bar),
        count_digits_after_decimal(δ)])
    Hminus_dum = (F .<= round(F_bar - δ; digits = sigdigs)) .& (F .> round(F_bar - 2*δ; digits = sigdigs))
    Hplus_dum = (F .>= round(F_bar + δ; digits = sigdigs)) .& (F .< round(F_bar + 2*δ; digits = sigdigs))
    B_dum = (F .> round(F_bar - δ; digits = sigdigs)) .& (F .< round(F_bar + δ; digits = sigdigs))

    #Get sums
    Hminus = sum(Hminus_dum)
    Hplus = sum(Hplus_dum)
    B_pre = sum(B_dum)

    #Get density heights by dividing by d
    Hminus_height = Hminus / δ
    Hplus_height = Hplus / δ

    #GENERATES B AND VARIABLE X APPEARING IN QUADRATIC
    B = B_pre - Hminus - Hplus
    #X = (F_bar*(Hplus_height - Hminus_height) / 2.0) - B
    #Solve for X using Quadratic Formula
    #θ = (-X + (X^2 + F_bar * Hminus_height * Hplus_height)^0.5) / (2.0 * F_bar * Hminus_height)

    a = F_bar * Hminus_height / 2.0
    b = F_bar * Hplus_height / 2.0 - F_bar * Hminus_height / 2.0 - B
    c = -F_bar * Hplus_height / 2.0

    θ = (-b + sqrt(b^2 - 4.0*a*c)) / (2.0 * a)

    println(Hplus_height)
    D = (Hplus_height) * ((1.0 + θ) / θ) * ((θ - 1.0) * F_bar)^2 / 2.0

    return θ, D

end


function bootstrap_θ(F::Vector{Float64}, d::Float64, F_bar::Float64, N_boot::Int64, ci_lvl::Float64)
    bs = bootstrap(X -> est_θ(X, d, F_bar), F, BasicSampling(N_boot))
    bci = confint(bs, BCaConfInt(ci_lvl))

    θ = bs.t0[1]
    θ_σ = stderror(bs)[1]
    θ_bias = bias(bs)[1]
    θ_bci = bci[1]
    θ_bs = bs.t1[1]

    return θ, θ_σ, θ_bias, θ_bci, θ_bs
end

function get_osr_limit(zone::AbstractString,s)
    if s < 1
        error("The number of stories, s, must be greater than or qual to 1.")
    end

    if zone == "R6"
        OSR_bar = 27.5 + (s-1)*0.5
    elseif zone in ["R7-1","R7-2","R7-3"]
        OSR_bar = 15.5 + (s-1)*0.5
    elseif zone[1:2] == "R8"
        OSR_bar = 5.9 + (s-1)*0.3
    elseif zone[1:2] == "R9"
        OSR_bar = 1.0 + (s-1)*0.4
    else
        error("Zoning District does not use Height Factor Regulations")
    end

    return OSR_bar / 100.0

end

function get_far_bar_hr(zone::AbstractString)
    if zone == "R6"
        FAR_bar = 2.43
    elseif zone in ["R7-1","R7-2","R7-3"]
        FAR_bar = 3.44
    elseif zone[1:2] == "R8"
        FAR_bar = 6.02
    elseif zone[1:2] == "R9"
        FAR_bar = 7.52
    else
        error("Zoning District does not use Height Factor Regulations")
    end
end

function is_hr_zoning(zone::AbstractString)

    if zone in ["R6","R7-1","R7-2","R7-3","R8","R9"]
        return true
    else
        return false
    end
end

function compute_max_bldg_front(zone::AbstractString, LotFront::Float64, s)
    osr_bar = get_osr_limit(zone,s)

    c1 = LotFront / (s * osr_bar + 1)
    far_bar = get_far_bar_hr(zone)
    c2 = far_bar * LotFront / s
    if c1 < c2
        return c1
    else
        return c2
    end
end

function compute_max_far_hr(zone::AbstractString, s)
    if s>0
        osr_bar = get_osr_limit(zone,s)
        far_bar = get_far_bar_hr(zone)
        return minimum([far_bar, s / (1 + s * osr_bar)])
    else 
        return missing
    end
end

function compute_s_star(zone::AbstractString, s_bar::Int64)

   return findmax([compute_max_far_hr(zone, s) for s = 1:s_bar])[2]

end

function compute_max_bldg_area(zone::AbstractString, LotDepth::Float64,
                                LotFront::Float64; s_bar = 100)

    s_star = compute_s_star(zone, s_bar)
    osr_limit = get_osr_limit(zone, s_star)
    far_limit = compute_max_far_hr(zone, s_star)
    Bf_star = compute_max_bldg_front(zone, LotFront, s_star)

    return Bf_star * LotDepth

end



function process_far_data(far_val::Float64)
    df_out = df[(df.MaxAllwFAR .== far_val) .& (df.YQ .== 2017.75), :]

    select!(df_out, Not([:YQ, :CondoNo, :CondoDum, :NumFloors_prev, :FAR_prev,
    :YearBuilt_prev, :YearAlter1_prev, :YearAlter2_prev, :NumFloors_chg,
    :FAR_chg, :YearBuilt_chg, :YearAlter1_chg, :YearAlter2_chg,
    :NumFloors_chg_sum, :NumFloors_chg_count, :ResidFAR, :CommFAR, :FacilFAR]))

    
    @rename! df_out :ZoneDist = :ZoneDist1
    @rename! df_out :Overlay = :Overlay1    
    @rename! df_out :FAR_bar = :MaxAllwFAR
    @rename! df_out :FAR_ratio = :FARRatio

    df_merge = df_fc[:,[:BoroCode, :Block, :Lot, :YQ, :YearBuilt, :YearAlter1, 
    :YearAlter2, :FAR_prev, :YearBuilt_prev, :YearAlter1_prev,
    :YearAlter2_prev]]
    @rename! df_merge :YQ_redev = :YQ

    function identify_nb(YearBuilt, YearBuilt_prev)
        if ismissing(YearBuilt_prev)
            return 1
        else
            if YearBuilt == YearBuilt_prev
                return 0
            else
                return 1
            end
        end
    end
    transform!(df_merge, [:YearBuilt, :YearBuilt_prev] => ByRow(identify_nb)
        => :NewBldg)
        
    
    df_merge = df_merge[:, Not([:YearBuilt, :YearAlter1, :YearAlter2,
        :YearAlter1_prev, :YearAlter2_prev])]
    
    df_out = leftjoin(df_out,df_merge, on = [:BoroCode, :Block, :Lot])

    return df_out

end

function process_zoning_data(zone::AbstractString, df::DataFrame)
    
    df_out = df[(df.ZoneDistYearBuilt .== zone) .& (df.Year .== 2017), :]

    select!(df_out, Not([:YQ, :CondoNo, :CondoDum, :ResidFAR, :CommFAR, :FacilFAR]))

    
    @rename! df_out :ZoneDist = :ZoneDist1
    @rename! df_out :Overlay = :Overlay1
    #@rename! df_out :FAR_bar = :MaxAllwFAR
    #@rename! df_out :FAR_ratio = :FARRatio

    #df_merge = df_fc[:,[:BoroCode, :Block, :Lot, :YQ, :YearBuilt, :YearAlter1, 
    #:YearAlter2, :FAR_prev, :YearBuilt_prev, :YearAlter1_prev,
    #:YearAlter2_prev]]
    #@rename! df_merge :YQ_redev = :YQ

    #function identify_nb(YearBuilt, YearBuilt_prev)
    #    if ismissing(YearBuilt_prev)
    #        return 1
    #    else
    #        if YearBuilt == YearBuilt_prev
    #            return 0
    #        else
    #            return 1
    #        end
    #    end
    #end
    #transform!(df_merge, [:YearBuilt, :YearBuilt_prev] => ByRow(identify_nb)
    #    => :NewBldg)
        
    
    #df_merge = df_merge[:, Not([:YearBuilt, :YearAlter1, :YearAlter2,
    #    :YearAlter1_prev, :YearAlter2_prev])]
    
    #df_out = leftjoin(df_out,df_merge, on = [:BoroCode, :Block, :Lot])

    df_out = df_out[.!ismissing.(df_out.YearBuilt), :]
    df_out = df_out[df_out.YearBuilt .>= 2000.0, :]

    @select!(df_out, :BoroCode, :Block, :Lot, :FAR, :NumFloors, :LotArea, :XCoord, :YCoord)
    return df_out

end

function generate_subsample(zonelist::Vector{String}, df::DataFrame)
    df_out = DataFrame()
    for (n,zonestr) in enumerate(zonelist)
        if n>1
            df_out = vcat(df_out, process_zoning_data(zonestr,df))
        else
            df_out = process_zoning_data(zonestr, df)
        end
    end

    return df_out
end

function generate_fbar_plus_allowance_res(zone_str)

    list_0p5 = ["R1-1","R1-2","R1-2A","R2", "R2A"]
    list_0p5_aa = ["R3-1","R3-2","R3A","R3X"]
    list_0p75 = ["R4","R4-1","R4A"]
    list_0p85 = ["R2X"]
    list_0p9 = ["R4B"]
    list_1p0 = ["M1-1"]
    list_1p1 = ["R5A"]
    list_1p135 = ["R5B"]
    list_1p25 = ["R5"]
    list_2p0 = ["R5D","M1-2","M1-4","M3-1","M3-2","M2-1","M2-3"]
    list_2p0_a = ["R6B"]
    list_3p0 = ["R7B"]
    list_3p0_a = ["R6A"]
    list_4p0 = ["R8B"]
    list_4p0_a = ["R7A"]
    list_4p2 = ["R7D"]
    list_5p0 = ["M1-3","M1-5","M2-2","M2-4"]
    list_5p0_a = ["R7X"]
    list_6p02 = ["R8A","R8X"]
    list_7p52 = ["R9A"]
    list_9p0_1 = ["R9D"]
    list_9p0_2 = ["R9X"]
    list_10p0 = ["M1-6"]
    list_10p0_a = ["R10A","R10X"]
    list_hf = ["R6","R7-1","R7-2","R8","R9","R10"]

    if zone_str in list_0p5
        return (0.5, 0.0, 0.0)
    elseif zone_str in list_0p5_aa
        return (0.5, 0.1, 0.0)
    elseif zone_str in list_0p75
        return (0.75, 0.15, 0.0)
    elseif zone_str in list_0p85
        return (0.85, 0.17, 0.0)
    elseif zone_str in list_0p9
        return (0.9, 0.0, 0.0)
    elseif zone_str in list_1p0
        return (1.0, 0.0, 0.0)
    elseif zone_str in list_1p1
        return (1.1, 0.0, 0.0)
    elseif zone_str in list_1p135
        return (1.135, 0.0, 0.0)
    elseif zone_str in list_1p25
        return (1.25, 0.0, 0.0)
    elseif zone_str in list_2p0
        return (2.0, 0.0, 0.0)
    elseif zone_str in list_2p0_a
        return (2.0, 0.2, 0.0)
    elseif zone_str in list_3p0
        return (3.0, 0.0, 0.0)
    elseif zone_str in list_3p0_a
        return (3.0, 0.0, 0.6)
    elseif zone_str in list_4p0
        return (4.0, 0.0, 0.0)
    elseif zone_str in list_4p0_a
        return (4.0, 0.0, 0.6)
    elseif zone_str in list_4p2
        return (4.2, 0.0, 1.4)
    elseif zone_str in list_5p0
        return (5.0, 0.0, 0.0)
    elseif zone_str in list_5p0_a
        return (5.0, 0.0, 1.0)
    elseif zone_str in list_6p02
        return (6.02, 0.0, 1.0)
    elseif zone_str in list_7p52
        return (7.52, 0.0, 0.98)
    elseif zone_str in list_9p0_1
        return (9.0, 0.0, 1.0)
    elseif zone_str in list_9p0_2
        return (9.0, 0.0, 0.7)
    elseif zone_str in list_10p0
        return (10.0, 0.0, 0.0)
    elseif zone_str in list_10p0_a
        return (10.0, 0.0, 2.0)
    elseif zone_str in list_hf
        return (0.0, 0.0, 0.0)
    else
        error("Zone string, $zone_str, not recognized")
    end
end

function generate_fbar_plus_allowance(zone_str)
    if string(zone_str[1]) != "C"
        res_out = generate_fbar_plus_allowance_res(zone_str)
        return res_out[1], res_out[2], res_out[3], 0.0
    else

        if zone_str == "C3"
            CFBar = 0.5
            res_str = "R3-2"
        elseif zone_str == "C3A"
            CFBar = 0.5
            res_str = "R3A"
        elseif zone_str == "C8-1"
            CFBar = 1.0
            return (0.0, 0.0, 0.0, CFBar)
        elseif zone_str == "C4-1"
            CFBar = 1.0
            res_str = "R5"
        elseif zone_str in ["C8-2","C8-3"]
            CFBar = 2.0
            return (0.0, 0.0, 0.0, CFBar)
        elseif zone_str in ["C1-6A","C2-6A"]
            CFBar = 2.0
            res_str = "R7A"
        elseif zone_str == "C1-7A"
            CFBar = 2.0
            res_str = "R8A"
        elseif zone_str in ["C2-7A", "C1-8A"]
            CFBar = 2.0
            res_str = "R9A"
        elseif zone_str in ["C1-8X","C2-7X"]
            CFBar = 2.0
            res_str = "R9X"
        elseif zone_str in ["C1-9A","C2-8A"]
            CFBar = 2.0
            res_str = "R10A"
        elseif zone_str in ["C4-2A","C4-3A"]
            CFBar = 3.0
            res_str = "R6A"
        elseif zone_str == "C4-D"
            CFBar = 3.4
            res_str = "R8A"
        elseif zone_str == "C4-6A"
            CFBar = 3.4
            res_str = "R10A"
        elseif zone_str in ["C4-4A","C4-5A"]
            CFBar = 4.0
            res_str = "R7A"
        elseif zone_str == "C4-5X"
            CFBar = 4.0
            res_str = "R7X"
        elseif zone_str == "C5-5"
            CFBar = 5.0
            res_str = "R10A"
        elseif zone_str == "C4-5D"
            CFBar = 4.2
            res_str = "R7D"
        elseif zone_str == "C8-4"
            CFBar = 5.0
            return (0.0, 0.0, 0.0, CFBar)
        elseif zone_str == "C6-2A"
            CFBar = 6.0
            res_str = "R8A"
        elseif zone_str == "C6-3A"
            CFBar = 6.0
            res_str = "R9A"
        elseif zone_str == "C6-3X"
            CFBar = 6.0
            res_str = "R9X"
        elseif zone_str == "C6-3D"
            CFBar = 9.0
            res_str = "R9D"
        elseif zone_str in ["C4-7A","C5-2A","C6-4A","C6-4X"]
            CFBar = 10.0
            res_str = "R10A"
        elseif zone_str == "C7"
            CFBar = 2.0
            return (0.0, 0.0, 0.0, CFBar)
        else
            error("Commercial zone $zone_str not found. Check for error.")
        end
        
        res_out = generate_fbar_plus_allowance_res(res_str)
        return res_out[1], res_out[2], res_out[3], CFBar
    end
end



function get_zoning_num(ZoneDist)
    
    match_obj1 = match(r"(?<=.)\d{1,2}", ZoneDist)

    if match_obj1 !== nothing
        return parse(Int64,match_obj1.match)
    else
        error("No match. Check erroneous zoning district.")
    end
end

function generate_fbar(zone_str)
    fbar,~,~,~ = generate_fbar_plus_allowance(zone_str)

    return fbar
end

function generate_attic_allowance(zone_str)
    ~,attic_allowance,~,~ = generate_fbar_plus_allowance(zone_str)

    return attic_allowance
end

function generate_ih_bonus(zone_str)
    ~,~,ih_bonus,~ = generate_fbar_plus_allowance(zone_str)

    return ih_bonus
end

function generate_CFbar(zone_str)
    ~,~,~,CFBar = generate_fbar_plus_allowance(zone_str)

    return CFBar
end



function generate_overlay_bonus(ZoneDist1, Overlay1)
    zone_use = string(ZoneDist1[1])
    zone_num = get_zoning_num(ZoneDist1)
    clist = ["C1-1","C1-2","C1-3","C1-4","C1-5","C2-1","C2-2","C2-3","C2-4","C2-5"]
    if ismissing(Overlay1) | (zone_use != "R")
        return 0.0
    end

    if Overlay1 in clist
        if zone_num < 6
            return 1.0
        else
            return 2.0
        end
    end

end

function generate_community_facility_fbar(zone_str)
    
    if string(zone_str[1]) != "R"
        return 0.0
    end

    if zone_str in ["R1-1","R1-2","R1-2A","R2","R2A","R2X","R3-1","R3-2","R3A","R3X"]
        return 1.0
    elseif zone_str in ["R4","R4-1","R4A","R4B","R5A","R5B","R5D","R6B","R5"]
        return 2.0
    elseif zone_str in ["R6A","R7B"]
        return 3.0
    elseif zone_str in ["R7A","R8B"]
        return 4.0
    elseif zone_str in ["R7D"]
        return 4.2
    elseif zone_str in ["R7X"]
        return 5.0
    elseif zone_str in ["R8X"]
        return 6.0
    elseif zone_str in ["R8A"]
        return 6.5
    elseif zone_str in ["R9A"]
        return 7.5
    elseif zone_str in ["R9D","R9X"]
        return 9.0
    elseif zone_str in ["R10X","R10A"]
        return 10.0
    elseif zone_str in ["R6","R7-1","R7-2","R8","R9","R10"]
        return 0.0
    else
        error("Zone string, $zone_str, not recognized")
    end
        
end

function drop_manufacturing_overlay(ZoneDist1, Overlay1)

    if startswith(ZoneDist1, "M") && !ismissing(Overlay1)
        return false
    else
        return true
    end
end

function drop_manufacturing_overlay(df::DataFrame)
    ZoneDist1 = df.ZoneDist1
    Overlay1 = df.Overlay1
    if (string(ZoneDist1[1])=="M") & !ismissing(Overlay1)
        return false
    else
        return true
    end
end

function filter_hf_commercial(zone_str)

    zone_list = ["C1-6","C1-7","C1-8","C1-9","C2-6","C2-7","C2-8","C4-2",
             "C4-3","C4-4","C4-5","C4-6","C4-7","C5-1","C5-2",
             "C5-3","C5-4","C5-5","C6-1","C6-2","C6-3","C6-4",
             "C6-5","C6-6","C6-7","C6-8","C6-9","C6-1A"]

    if string(zone_str[1]) == "C"
        if zone_str in zone_list
            return false
        else
            return true
        end
    else
        return true
    end
end

def_condodum = function(x)
    if ismissing(x)
        return 0
    elseif x[1] == 'R'
        return 1
    else
        return 0
    end
end

function correct_year_built(x)
    if (typeof(x) == String) 
        return parse(Int64, x)
    else
        return x
    end
end

function filter_bldg_class(x, bldg_class_filter)
    if string(unwrap(x)[1]) in bldg_class_filter
        x[1] in bldg_class_filter
        return false
    else
        return true
    end
end

function correct_max_allw_far(x)
    if typeof(x) == String
        return parse(Float64, x)
    else
        return x
    end
end

function xform_maxfar(max_far, com_area, res_area, resid_far, com_far)
    if ismissing(max_far)
        if com_area > res_area
            return com_far
        else
            return resid_far
        end
    else
        return max_far
    end
end

function correct_boro_code(x)
    if (typeof(x) == String)  .|| (typeof(x) == String15)
        return parse(Int64, x)
    else
        return x
    end
end


function find_zone_dist(BoroCode::Int64, Block::Int64, Lot::Int64, year::Int64,
    df::DataFrame)
    
    df_subset = @rsubset df begin
        :BoroCode == BoroCode
        :Block == Block
        :Lot == Lot
    end

    min_year = minimum(df_subset.Year)
    if year >= min_year
        tmp =  df_subset.ZoneDist1[df_subset.Year .== year]
        while isempty(tmp) && year <= 2017
            year += 1
            tmp = df_subset.ZoneDist1[df_subset.Year .== year]
        end
        return tmp[1]
    else
        return df_subset.ZoneDist1[1]
    end

end

function zoning_dist_year_built(BoroCode::Int64, Block::Int64, Lot::Int64, YearBuilt::Int64,
    df::DataFrame)

    return find_zone_dist(BoroCode, Block, Lot, maximum([YearBuilt, 2002]), df)

end

function is_mixed_zone(ZoneDist::String)
    return !occursin("/", ZoneDist)  # Return true if ZoneDist does NOT contain "/"
end

export clean_raw_pluto_data, transform_cleaned_pluto_data, generate_analysis_subsamples,
    create_sample_map, generate_tables_figures, est_θ, bootstrap_θ

end