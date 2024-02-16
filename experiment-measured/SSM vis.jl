using JLD2
using MAT
using Plots
default(label="", linewidth=3, margin=15Plots.pt)

## Load data
jldopen("results/")