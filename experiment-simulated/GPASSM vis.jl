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
@unpack true_mcp_1, true_mcp_2, true_mcp_3, true_τ_a, tsteps, Δt, fitx_m, fitx_v, states, measurements, inputs, tsteps_val, sim_states, states_val, μ1,Σ1, μ2,Σ2, μ3,Σ3, ρ1_m,ρ1_s, ρ2_m,ρ2_s, ρ3_m,ρ3_s, l_star, γ_star = file

# True nonlinear convection function
r(T_i,T_a) = (T_a - T_i)^3 ./ 100

## Visualize data set

hline([true_τ_a], color="black", linestyle=:dash, label=L"$T_a$")
scatter!(tsteps ./ 60, 
      transpose(measurements), 
      xlabel="time [min]",
      ylabel="temperature [C]",
      xticks = [0, 5, 10, 15],
      ylims = (19,37),
      yticks = [20,25,30,35],
      grid=true,
      alpha = 1.,
      markersize=2,
      markerstrokecolors = ["red" "blue" "orange"],
      markercolors = ["red" "blue" "orange"], 
      labels = [L"$T_1$" L"$T_2$" L"$T_3$"],
      legendfontsize=10,
      title="simulated data",
      size=(400,250),
)

savefig("experiment-simulated/figures/simulatedexp-temperatures.png")

plot(tsteps ./ 60, 
      transpose(inputs), 
      xlabel="time [min]",
      ylabel="input [W]",
      xticks = [0, 5, 10, 15],
      grid=true,
      linecolors = ["darkred" "darkblue" "darkorange"],
      linestyles = :dot,
      labels = [L"$u_1$" L"$u_2$" L"$u_3$"],
      legendfontsize=10,
      legend=:right,
      size=(400,150),
)

savefig("experiment-simulated/figures/simulatedexp-inputs.png")


## Visualize nonlinearity estimates

plot(xlabel="temperature [C]",
     ylabel="r(T₁,Tₐ) [C]",
     size=(400,300),
     legend=:bottomleft,
     legendfontsize=10,
     )
plot!(states[1,:], r.(states[1,:],true_τ_a), linewidth=5, linestyle=:solid, color="darkred", alpha=1., label="true")
scatter!(fitx_m[1,:], fitx_m[4,:], alpha=0.5, markerstrokewidth=0, markersize=fitx_v[4,:]/2, color="red", label="GP states")
plot!(sort(fitx_m[1,:]), ρ1_m, ribbon=ρ1_s, color="black", linestyle=:dash, fillalpha=0.2, linewidth=4, label="polynomial fit")

savefig("experiment-simulated/figures/VL+GPASSM-block1_fnest.png")

plot(xlabel="temperature [C]",
     ylabel="r(T₂,Tₐ) [C]",
     legend=:bottomleft,
     legendfontsize=10,
     size=(400,300),
     )
plot!(states[2,:], r.(states[2,:],true_τ_a), linewidth=5, color="lightblue", linestyle=:solid, alpha=1., label="true")
scatter!(fitx_m[2,:], fitx_m[5,:], alpha=0.5, markerstrokewidth=0, markersize=fitx_v[5,:]/2, color="blue", label="GP states")
plot!(sort(fitx_m[2,:]), ρ2_m, ribbon=ρ2_s, color="black", linestyle=:dash, fillalpha=0.2, linewidth=4, label="polynomial fit")

savefig("experiment-simulated/figures/VL+GPASSM-block2_fnest.png")

plot(xlabel="temperature [C]",
     ylabel="r(T₃,Tₐ) [C]",
     legend=:bottomleft,
     legendfontsize=10,
     size=(400,300),
     )
plot!(states[3,:], r.(states[3,:],true_τ_a), linewidth=5, color="darkorange", linestyle=:solid, alpha=1., label="true")
scatter!(fitx_m[3,:], fitx_m[6,:], alpha=0.5, markerstrokewidth=0, markersize=fitx_v[6,:] ./2, color="orange", label="GP states")
plot!(sort(fitx_m[3,:]), ρ3_m, ribbon=ρ3_s, color="black", linestyle=:dash, fillalpha=0.2, linewidth=4, label="polynomial fit")

savefig("experiment-simulated/figures/VL+GPASSM-block3_fnest.png")

## Visualize simulation error

SMSE_GPSSM = mean((sim_states - states_val).^2)

scatter(tsteps_val[1:10:end] ./ 60, 
     transpose(states_val[:,1:10:end]), 
     markersize=5,
     markerstrokewidth=1,
     markercolors = ["red" "blue" "orange"], 
     marker=".",
)
scatter!([1e8], [1e8], ylims=[22, 35], xlims=[0, 16], marker='.', color="gray", label="system")
plot!(tsteps_val ./ 60,
     sim_states',
     # ribbon=,
     linewidth=4,
     linecolors = ["red" "blue" "orange"], 
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
     xlims=[0, 16],
     size=(400,300),
     legend=:topright,
)

savefig("experiment-simulated/figures/LGPASSM-simulation-errors.png")