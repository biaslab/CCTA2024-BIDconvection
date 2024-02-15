using JLD2
using MAT
using Plots
default(label="", linewidth=3, margin=15Plots.pt)

## Load data
jldopen("results/")

## Visualize white-box estimates

plot(tsteps ./ 60,
     fitx_m';
     ribbon=sqrt.(fitx_v)',
     legend = true, 
     # title="GPSSM Smoothed MSE = $FMSE_GPSSM",
     linecolors = ["red" "blue" "orange"], 
     fillcolors = ["red" "blue" "orange"], 
     labels = [L"$τ_1$ estimated" L"$τ_2$ estimated" L"$τ_3$ estimated"],
     xlabel = "time [min]", 
     ylabel = "temperature [C]",
     size=(500,400)
)
# plot!(tsteps, 
#       transpose(states), 
#       alpha = 0.3,
#       linecolors = ["red" "blue" "orange"], 
#       linestyle = :dash,
#       labels = [L"$τ_1$ true" L"$τ_2$ true" L"$τ_3$ true"],
# )
# plot!(tsteps, states', linecolors = ["red" "blue" "orange"], labels = [L"$τ_1$ true" L"$τ_2$ true" L"$τ_3$ true"], )
scatter!(tsteps ./ 60, measurements', markercolors = ["red" "blue" "orange"], labels = [L"$τ_1$ observed" L"$τ_2$ observed" L"$τ_3$ observed"])
# hline!([true_τ_a], color="black", linewidth=1, linestyle=:dot, label="ambient temp τ_a")

savefig("figures/SSM+states.png")

# Post-measurement residuals
residuals = measurements - C*fitx_m

plot(tsteps, residuals[1,:], color="red")
plot!(tsteps, residuals[2,:], color="blue")
plot!(tsteps, residuals[3,:], color="orange")

savefig("figures/SSM+GPr-residuals.png")


## Visualize nonlinearity estimates

plot(xlabel="",
     ylabel="convection [C]",
     # yticks=(range(-1000,stop=0.0, length=5), round.(range(-1000,stop=0.0, length=5).*Δt./true_mcp_1, digits=2)),
     size=(500,250),
     legend=:bottomleft,
     )
plot!(states[1,:], r.(states[1,:],1,true_τ_a) ./ true_mcp_1, linewidth=5, color="black", alpha=1., label="true")
# plot!(fitx_m[1,:], fitx_m[4,:], ribbon=sqrt.(fitx_v[4,:]), color="red", alpha=0.5, fillalpha=0.5, label="GP")
scatter!(fitx_m[1,:], rGP_1, alpha=1., markerstrokewidth=0, markersize=2, label="GP fit", color="red")
# scatter!(fitx_m[1,:], samples_f1[2:end,:]', alpha=0.05, markerstrokewidth=0, markersize=3, label="", color="red")
plot!(sort(fitx_m[1,:]), x -> f1_pol3(x), linewidth=4, color="darkred", linestyle=:dash, label="polynomial")

savefig("figures/SSM+GPr-block1_fnest.png")

plot(xlabel="",
     ylabel="convection [C]",
     # yticks=(range(-1000,stop=0.0, length=5), round.(range(-1000,stop=0.0, length=5).*Δt./true_mcp_1, digits=2)),
     size=(500,200),
     legend=:false,
     )
plot!(states[2,:], r.(states[2,:],1,true_τ_a) ./ true_mcp_2, linewidth=5, color="black", alpha=1., label="true")
# plot!(fitx_m[1,:], fitx_m[4,:], ribbon=sqrt.(fitx_v[4,:]), color="red", alpha=0.5, fillalpha=0.5, label="GP")
scatter!(fitx_m[2,:], rGP_2, alpha=1., markerstrokewidth=0, markersize=2, label="GP fit", color="blue")
# scatter!(fitx_m[2,:], samples_f2[2:end,:]', alpha=0.05, markerstrokewidth=0, markersize=3, label="", color="blue")
plot!(sort(fitx_m[2,:]), x -> f2_pol3(x), linewidth=4, color="darkblue", linestyle=:dash, label="polynomial")

savefig("figures/SSM+GPr-block2_fnest.png")

plot(xlabel="temperature [C]",
     ylabel="convection [C]",
     # yticks=(range(-1000,stop=0.0, length=5), round.(range(-1000,stop=0.0, length=5).*Δt./true_mcp_1, digits=2)),
     size=(500,200),
     legend=:false,
     )
plot!(states[3,:], r.(states[3,:],1,true_τ_a) ./ true_mcp_3, linewidth=5, color="black", alpha=1., label="true")
# plot!(fitx_m[1,:], fitx_m[4,:], ribbon=sqrt.(fitx_v[4,:]), color="red", alpha=0.5, fillalpha=0.5, label="GP")
scatter!(fitx_m[3,:], rGP_3, alpha=1., markerstrokewidth=0, markersize=2, label="GP fit", color="orange")
# scatter!(fitx_m[2,:], samples_f2[2:end,:]', alpha=0.05, markerstrokewidth=0, markersize=3, label="", color="blue")
plot!(sort(fitx_m[3,:]), x -> f3_pol3(x), linewidth=4, color="darkorange", linestyle=:dash, label="polynomial")

savefig("figures/SSM+GPr-block3_fnest.png")

## Visualize simulation error

SMSE_GPSSM = mean((sim_states - states_val).^2)

plot(tsteps_val ./ 60,
     sim_states';
     legend = true, 
     linewidth=4,
    #  title = "GPSSM Simulation MSE = $SMSE_GPSSM",
     linecolors = ["red" "blue" "orange"], 
     fillcolors = ["red" "blue" "orange"], 
     labels = [L"$T_1$ identified" L"$T_2$ identified" L"$T_3$ identified"],
     xlabel = "time [min]", 
     ylabel = "temperature [C]",
     size=(500,250)
)
plot!(tsteps_val ./ 60, 
      transpose(states_val), 
      linewidth=4,
      alpha = 1.,
      linecolors = ["red" "blue" "orange"], 
      linestyle = :dash,
      labels = [L"$T_1$ system" L"$T_2$ system" L"$T_3$ system"],
)

savefig("figures/NONLCONV_SSM+GPr-simulations.png")