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
@unpack M,F,A,B,C,Q,R, measurements_ss, inputs_ss, tsteps, Δt, fitx_m, fitx_v, h_pol, λ_star, γ_star, sim_temps = file

## Plot measured temperatures

vis_ss = 1
vis_ix = [3,9,12]

lcolors = ["red" "blue" "orange"]
labelsm = hcat([L"$T_{%$k}$" for k in vis_ix]...)

plt1 = plot(tsteps[2:vis_ss:end] ./ 3600, 
            measurements_ss[2:vis_ss:end,vis_ix],
            linecolors = lcolors, 
            linewidth=4,
            labels = labelsm,
            xscale=:log2,
            ylabel="temperature [C]",
            legend=:topleft,
            # xticks=([log10(8.), log10(16.), log10(24.)], [8., 16., 24.]),
            ylims=[20., 50.],
            grid=true)
plt2 = plot(tsteps[2:vis_ss:end] ./ 3600, 
            inputs_ss[2:vis_ss:end,1],
            xlabel="time [h]",
            color="black",
            yticks=[23, 24],
            ylims=[23, 24],
            label=L"$T_a$",
            xscale=:log10,
            ylabel="temperature [C]",
            legend=:top)            
plot(plt1,plt2, layout=grid(2,1, heights=[0.6, 0.4]), margin=0.1Plots.mm, size=(500,350))

savefig("figures/measured-temperatures.png")

## Visualize as a function of both T_i and T_a

ntc_ix = 3
vis_ss = 10

scatter3d(inputs_ss[2:vis_ss:end,1],
         fitx_m[ntc_ix,2:vis_ss:end],
         fitx_m[12+ntc_ix,2:vis_ss:end],
         legend = true, 
        color = "red",
#      labels = "GP NTC#$ntc_ix",
        title ="NTC#$ntc_ix",
        ylabel = "Block temperature [C]",
        xlabel = "Ambient temperature [C]",
        zlabel = "Change in temperature [ΔC]",
      #  camera = (90,20),
        size=(500,500)
)

## Visualize GP fits as a function of T_i

ntc_ix = 1

plot(xlabel="temperature [C]",
     ylabel="change in temperature [ΔC]",
     title ="NTC#$ntc_ix"
)
scatter!(fitx_m[ntc_ix,:], fitx_m[N+ntc_ix,:], markersize=1e-2*sqrt.(fitx_v[N+ntc_ix,:]*Δt./M[ntc_ix,ntc_ix]), alpha=0.1, color="red")
plot!(sort(fitx_m[ntc_ix,:]), x -> h_pol[ntc_ix](x), color="black", linestyle=:dash, label="polynomial fit")

savefig("figures/GP-NTC$ntc_ix-fn.png")