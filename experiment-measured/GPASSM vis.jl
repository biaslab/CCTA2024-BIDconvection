using StatsBase
using LinearAlgebra
using JLD2
using UnPack
using Polynomials
using LaTeXStrings
using MAT
using Plots
default(label="", linewidth=3, margin=5Plots.pt)

## Load data
file = jldopen("experiment-measured/results/GPASSM.jld2")
@unpack M,F,A,B,C,Q,R, measurements_ss, inputs_ss, N, tsteps, Δt, fitx_m, fitx_v, h_pol, λ_star, γ_star, sim_temps = file

## Plot measured temperatures

vis_ss = 1
vis_ix = [3,9,12]

lcolors = ["red" "blue" "orange"]
labelsm = hcat([L"$T_{%$k}$" for k in vis_ix]...)

plt1 = plot(tsteps ./ 3600, 
            measurements_ss[:,vis_ix],
            linecolors = lcolors, 
            linewidth=4,
            labels = labelsm,
        #     xscale=:log2,
            ylabel="temp. [C]",
            legend=:bottomright,
            # xticks=([log10(8.), log10(16.), log10(24.)], [8., 16., 24.]),
            ylims=[20., 50.],
            guidefontsize=8,
            legendfontsize=8,
            grid=true)
plt2 = plot(tsteps ./ 3600, 
            inputs_ss[:,1],
            xlabel="time [h]",
            color="black",
            yticks=[23, 24],
            ylims=[23, 24],
            label=L"$T_a$",
        #     xscale=:log10,
            ylabel="temp.[C]",
            guidefontsize=8,
            legendfontsize=8,
            legend=:topright)            
plot(plt1,plt2, layout=grid(2,1, heights=[0.6, 0.4]), margin=1.5Plots.mm, size=(500,250))

savefig("experiment-measured/figures/measured-temperatures.png")

## Visualize as a function of both T_i and T_a

ntc_ix = 3
vis_ss = 10

scatter3d(inputs_ss[2:vis_ss:end,1],
          fitx_m[ntc_ix,2:vis_ss:end],
          fitx_m[12+ntc_ix,2:vis_ss:end],
          color = "red",
          markersize=3,
          markerstrokecolor="red",
          alpha=0.1,
#        labels = "GP NTC#$ntc_ix",
          title ="NTC#$ntc_ix",
          ylabel = "T₃ [C]",
          yrotation = 30,
          xlabel = "Tₐ [C]",
          xrotation = -25,
          zlabel = "temperature change [ΔC]",
          camera = (60,20),
          size=(600,400),
          xlims=[23.2, 24.2],
          xticks=[23.5,24],
          ylims=[21, 46],
          yticks=[25, 30, 35, 40, 45],
)

savefig("experiment-measured/figures/NTC#$ntc_ix-GP-Ti+Ta.png")

## Visualize GP fits as a function of T_i

ntc_ix = 3
mcolor = "red"
plot(xlabel="temperature [C]",
     ylabel=L"change [ΔC $\times 10^{-5}$]",
#      title ="NTC#$ntc_ix",
#      yticks=[-6e-5, -4e-5, -2e-5, 0.0, 2e-5],
#      yformatter=:scientific,
     size=(500,200),
)
scatter!([fitx_m[ntc_ix,1]], [fitx_m[N+ntc_ix,1]].*1e5, markersize=10_000fitx_v[N+ntc_ix,1], alpha=1.0, markerstrokecolor=mcolor, color=mcolor, label="GP states")
scatter!(fitx_m[ntc_ix,:], fitx_m[N+ntc_ix,:].*1e5, markersize=20_000fitx_v[N+ntc_ix,:], alpha=0.1, markerstrokecolor=mcolor, color=mcolor, label="")
plot!(sort(fitx_m[ntc_ix,:]), x -> h_pol[ntc_ix](x)*1e5, color="black", linestyle=:dash, label="polynomial fit")

savefig("experiment-measured/figures/GP-NTC$ntc_ix-fn.png")

ntc_ix = 9
mcolor = "blue"
plot(xlabel="temperature [C]",
     ylabel=L"change [ΔC $\times 10^{-4}$]",
#      title ="NTC#$ntc_ix",
#      yticks=[-6e-5, -4e-5, -2e-5, 0.0, 2e-5],
#      yformatter=:scientific,
     size=(500,200),
)
scatter!([fitx_m[ntc_ix,1]], [fitx_m[N+ntc_ix,1]]*1e4, markersize=10_000fitx_v[N+ntc_ix,1], alpha=1.0, markerstrokecolor=mcolor, color=mcolor, label="GP states")
scatter!(fitx_m[ntc_ix,:], fitx_m[N+ntc_ix,:]*1e4, markersize=20_000fitx_v[N+ntc_ix,:], alpha=0.1, markerstrokecolor=mcolor, color=mcolor, label="")
plot!(sort(fitx_m[ntc_ix,:]), x -> h_pol[ntc_ix](x)*1e4, color="black", linestyle=:dash, label="polynomial fit")

savefig("experiment-measured/figures/GP-NTC$ntc_ix-fn.png")

ntc_ix = 12
mcolor = "orange"
plot(xlabel="temperature [C]",
     ylabel=L"change [ΔC $\times 10^{-4}$]",
#      title ="NTC#$ntc_ix",
#      yticks=[-6e-5, -4e-5, -2e-5, 0.0, 2e-5],
#      yformatter=:scientific,
     size=(500,200),
)
scatter!([fitx_m[ntc_ix,1]], [fitx_m[N+ntc_ix,1]]*1e4, markersize=10_000fitx_v[N+ntc_ix,1], alpha=1.0, markerstrokecolor=mcolor, color=mcolor, label="GP states")
scatter!(fitx_m[ntc_ix,:], fitx_m[N+ntc_ix,:]*1e4, markersize=20_000fitx_v[N+ntc_ix,:], alpha=0.1, markerstrokecolor=mcolor, color=mcolor, label="")
plot!(sort(fitx_m[ntc_ix,:]), x -> h_pol[ntc_ix](x)*1e4, color="black", linestyle=:dash, label="polynomial fit")

savefig("experiment-measured/figures/GP-NTC$ntc_ix-fn.png")

## Visualize simulation errors

ntc_ix = 3
plot(tsteps ./ 3600,
     measurements_ss[:,ntc_ix] .- sim_temps[:,ntc_ix],
     color="red", 
     # labels = labelss,
     xlabel = "time [h]", 
     ylabel = "error [C]",
     label = "GPLFM",
     size =(500,200),
#      title="MSE = $SIM_MSE",
     # xscale=:log10,
)

savefig("experiment-measured/figures/NTC$ntc_ix-simerror.png")

ntc_ix = 9
plot(tsteps ./ 3600,
     measurements_ss[:,ntc_ix] .- sim_temps[:,ntc_ix],
     color="blue", 
     # labels = labelss,
     xlabel = "time [h]", 
     ylabel = "error [C]",
     size=(500,200),
     label = "GPLFM",
#      title="MSE = $SIM_MSE",
     # xscale=:log10,
)

savefig("experiment-measured/figures/NTC$ntc_ix-simerror.png")

ntc_ix = 12
plot(tsteps ./ 3600,
     measurements_ss[:,ntc_ix] .- sim_temps[:,ntc_ix],
     color="orange", 
     # labels = labelss,
     xlabel = "time [h]", 
     ylabel = "error [C]",
     size=(500,200),
     label="GPLFM",
#      title="MSE = $SIM_MSE",
     # xscale=:log10,
)

savefig("experiment-measured/figures/NTC$ntc_ix-simerror.png")