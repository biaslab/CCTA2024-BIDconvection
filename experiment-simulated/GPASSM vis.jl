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
file = load("experiment-simulated/results/LGPASSM.jld2")
@unpack true_mcp_1, true_mcp_2, true_mcp_3, true_τ_a, tsteps, Δt, fitx_m, fitx_v, states, measurements, inputs, tsteps_val, sim_states, states_val, f1_pol3, f2_pol3, f3_pol3, l_star, γ_star = file
f1 = Polynomial(f1_pol3.coeffs)
f2 = Polynomial(f2_pol3.coeffs)
f3 = Polynomial(f3_pol3.coeffs)

r(T_i,T_a) = (T_a - T_i)^3 ./ 100

## Visualize data set

scatter(tsteps ./ 60, 
      transpose(measurements), 
      xlabel="time [min]",
      ylabel="temperature [C]",
      xticks = [0, 5, 10, 15],
      grid=true,
      alpha = 1.,
      markersize=2,
      markerstrokecolors = ["red" "blue" "orange"],
      markercolors = ["red" "blue" "orange"], 
      labels = [L"$T_1$" L"$T_2$" L"$T_3$"],
      legendfontsize=10,
      size=(400,300),
)

savefig("experiment-simulated/figures/simulated-temperatures.png")


## Visualize nonlinearity estimates

plot(xlabel="temperature [C]",
     ylabel="temperature change [ΔC]",
     # yticks=(range(-1000,stop=0.0, length=5), round.(range(-1000,stop=0.0, length=5).*Δt./true_mcp_1, digits=2)),
     # ylims=(-1000,0)
     size=(400,300),
     legend=:bottomleft,
     legendfontsize=10,
     )
plot!(states[1,:], r.(states[1,:],true_τ_a), linewidth=5, color="black", alpha=0.5, label="true r(T₁,Tₐ)")
# scatter!(fitx_m[1,:], samples_f1', alpha=0.05, markerstrokewidth=0, markersize=3, color="red")
scatter!([fitx_m[1,1]], [fitx_m[4,1]], alpha=1., markerstrokewidth=0, markersize=0.001, color="red", label="GP states")
scatter!(fitx_m[1,:], fitx_m[4,:], alpha=0.1, markerstrokewidth=0, markersize=sqrt.(fitx_v[4,:]) ./10, color="red", label="")
plot!(sort(fitx_m[1,:]), x -> f1(x), color="black", linewidth=4, linestyle=:dash, label="polynomial fit")

savefig("experiment-simulated/figures/VL+GPASSM-block1_fnest.png")

plot(xlabel="temperature [C]",
     ylabel="temperature change [ΔC]",
     # yticks=(range(-250,stop=50., length=4), round.(range(-250,stop=50., length=4).*Δt./true_mcp_2, digits=2)),
     # ylims=(-250,50),
     legend=:bottomleft,
     legendfontsize=10,
     size=(400,300),
     )
plot!(states[2,:], r.(states[2,:],true_τ_a), linewidth=5,color="black", alpha=0.5, label="true r(T₂,Tₐ)")
# scatter!(fitx_m[2,:], samples_f2', alpha=0.05, markerstrokewidth=0, markersize=3, color="blue")
scatter!([fitx_m[2,1]], [fitx_m[5,1]], alpha=1., markerstrokewidth=0, markersize=0.001, color="blue", label="GP states")
scatter!(fitx_m[2,:], fitx_m[5,:], alpha=0.1, markerstrokewidth=0, markersize=sqrt.(fitx_v[5,:]) ./10, color="blue")
plot!(fitx_m[2,:], x -> f2(x), color="black", linewidth=4, linestyle=:dash, label="polynomial fit")

savefig("experiment-simulated/figures/VL+GPASSM-block2_fnest.png")

plot(xlabel="temperature [C]",
     ylabel="temperature change [ΔC]",
     # ylims=(-250,50),
     legend=:bottomleft,
     legendfontsize=10,
     size=(400,300),
     )
plot!(states[3,:], r.(states[3,:],true_τ_a), linewidth=5, color="black", alpha=0.5, label="true r(T₃,Tₐ)")
# scatter!(fitx_m[3,:], samples_f3', alpha=0.05, markerstrokewidth=0, markersize=3, color="orange")
scatter!([fitx_m[3,1]], [fitx_m[6,1]], alpha=1., markerstrokewidth=0, markersize=0.001, color="orange", label="GP states")
scatter!(fitx_m[3,:], fitx_m[6,:], alpha=0.1, markerstrokewidth=0, markersize=sqrt.(fitx_v[6,:]) ./10, color="orange")
plot!(fitx_m[3,:], x -> f3(x), color="black", linewidth=4, linestyle=:dash, label="polynomial fit")

savefig("experiment-simulated/figures/VL+GPASSM-block3_fnest.png")

## Visualize simulation error

SMSE_GPSSM = mean((sim_states - states_val).^2)

scatter(tsteps_val[1:10:end] ./ 60, 
         transpose(states_val[:,1:10:end]), 
         markersize=10,
         markercolors = ["red" "blue" "orange"], 
         marker=:rect,
     #  labels = [L"$T_1$ true" L"$T_2$ true" L"$T_3$ true"],
)
scatter!([1e8], [1e8], ylims=[22, 35], xlims=[0, 15], marker=:x, color="gray", label="system")
plot!(tsteps_val ./ 60,
        sim_states';
        linewidth=4,
        linecolors = ["red" "blue" "orange"], 
        linestroke
     # labels = ["simulated T₁" "simulated T₂" "simulatd T₃"],
        xlabel = "time [min]", 
        ylabel = "temperature [C]",
        size=(400,300),
        legend=:topleft,
)
plot!([1e8], [1e8], color="black", linestyle=:solid, label="identified")

savefig("experiment-simulated/figures/LGPASSM-simulations.png")

## Simulation error
plot(tsteps_val ./ 60,
      transpose(states_val .- sim_states);
      linewidth=4,
      linecolors = ["red" "blue" "orange"], 
     labels = [L"T_1" L"T_2" L"T_3"],
     xlabel = "time [min]", 
     ylabel = "error [C]",
     size=(400,300),
     legend=:topright,
     legendfontsize=10,
     guidefontsize=10,
)

savefig("experiment-simulated/figures/LGPASSM-simulation-errors.png")