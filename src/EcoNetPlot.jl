module EcoNetPlot

# using EcologicalNetworksDynamics
using EcologicalNetworksDynamics
using GraphMakie
using CairoMakie
using Graphs
using UMAP
using Statistics
using LinearAlgebra
using Serialization

const END = EcologicalNetworksDynamics

include("plot.jl")
include("chilean-web.jl")

end

