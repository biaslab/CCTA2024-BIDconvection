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
file = load("experiment-simulated/results/SSM+GPr.jld2") 
@unpack true_mcp_1, true_mcp_2, true_mcp_3, true_τ_a, tsteps, Δt, fitx_m, fitx_v, states, measurements, inputs, tsteps_val, sim_states, states_val, f1_pol3, f2_pol3, f3_pol3, λ_star, γ_star, residuals = file
f1 = Polynomial(f1_pol3.coeffs)
f2 = Polynomial(f2_pol3.coeffs)
f3 = Polynomial(f3_pol3.coeffs)

r(T_i,T_a) = (T_a - T_i)^3 ./ 100

## Visualize white-box estimates

plot(tsteps ./ 60,
     fitx_m';
     ribbon=sqrt.(fitx_v'),
     legend = true, 
     linecolors = ["red" "blue" "orange"], 
     fillcolors = ["red" "blue" "orange"], 
     labels = [L"$T_1$ estimates" L"$T_2$ estimates" L"$T_3$ estimates"],
     xlabel = "time [min]", 
     ylabel = "temperature [C]",
     size=(400,300)
)
scatter!([tsteps[1]], measurements[:,1]', alpha=1., markercolors = ["red" "blue" "orange"], markerstrokecolors = ["red" "blue" "orange"], labels = [L"$τ_1$ measured" L"$τ_2$ measured" L"$τ_3$ measured"])
scatter!(tsteps ./ 60, measurements', alpha=0.1, markercolors = ["red" "blue" "orange"], markerstrokecolors = ["red" "blue" "orange"], labels ="")

savefig("experiment-simulated/figures/SSM+states.png")

# Post-measurement residuals
residuals = measurements - fitx_m

plot(tsteps, residuals[1,:], color="red")
plot!(tsteps, residuals[2,:], color="blue")
plot!(tsteps, residuals[3,:], color="orange")

# savefig("experiment-simulated/figures/SSM+GPr-residuals.png")


## Visualize residual estimates

plot(xlabel="temperature [C]",
     ylabel="residuals",
     size=(400,300),
     legend=:bottomleft,
     )
scatter!(fitx_m[1,:], residuals[1,:], alpha=1., markerstrokewidth=0, markersize=4, label="SSM residuals", color="red")
plot!(sort(fitx_m[1,:]), x -> f1_pol3(x), linewidth=4, color="black", linestyle=:dash, label="polynomial fit")

savefig("experiment-simulated/figures/SSM+GPr-block1_fnest.png")

plot(xlabel="temperature [C]",
     ylabel="residuals",
     size=(400,300),
     legend=:bottomleft,
     )
scatter!(fitx_m[2,:], residuals[2,:], alpha=1., markerstrokewidth=0, markersize=4, label="SSM residuals", color="blue")
plot!(sort(fitx_m[2,:]), x -> f2_pol3(x), linewidth=4, color="black", linestyle=:dash, label="polynomial fit")

savefig("experiment-simulated/figures/SSM+GPr-block2_fnest.png")

plot(xlabel="temperature [C]",
     ylabel="residuals",
     size=(400,300),
     legend=:bottomleft,
     )
scatter!(fitx_m[3,:], residuals[3,:], alpha=1., markerstrokewidth=0, markersize=4, label="SSM residuals", color="orange")
plot!(sort(fitx_m[3,:]), x -> f3_pol3(x), linewidth=4, color="black", linestyle=:dash, label="polynomial fit")

savefig("experiment-simulated/figures/SSM+GPr-block3_fnest.png")

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
     linewidth=4,
     linecolors = ["red" "blue" "orange"], 
     xlabel = "time [min]", 
     ylabel = "temperature [C]",
     size=(400,300),
     legend=:topleft,
)
plot!([1e8], [1e8], color="black", linestyle=:solid, label="identified")

savefig("experiment-simulated/figures/SSM+resid-pol-simulations.png")

## Simulation errors

plot(tsteps_val ./ 60,
      transpose(states_val .- sim_states);
      linewidth=4,
      linecolors = ["red" "blue" "orange"], 
     labels = [L"T_1" L"T_2" L"T_3"],
     xlabel = "time [min]", 
     ylabel = "error [C]",
     xlims=[0, 16],
     size=(400,300),
     legend=:topleft,
)

savefig("experiment-simulated/figures/SSM+resid-pol-simulation-errors.png")