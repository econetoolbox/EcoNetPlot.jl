# Plotting ecological networks

A plotting package for ecological networks (food webs and more) powered by [`GraphMakie`](https://github.com/MakieOrg/GraphMakie.jl) 
and [`UMAP`](https://github.com/dillondaudert/UMAP.jl).
Direct integration with [`EcologicalNetworksDynamics`](https://github.com/econetoolbox/EcologicalNetworksDynamics.jl).

### Installation

Install with

```julia
using Pkg; Pkg.add(url = "https://github.com/econetoolbox/EcoNetPlot.jl")
```


### Plotting from an adjacency matrix

You have an adjacency matrix which represents an ecological networks, say, a food web.
Elements `A[i,j] = 1` means that species i
eats species j.
Adjacency matrix are useful tool to encode networks, but they are not very visual.
The present package allows you to plot networks from adjacency matrix.

For example, let's that we have two herbivores (2 and 3) eating the same plant (1).
To plot this network, we simple need to write the corresponding adjacency matrix
and give it to the function [`plot_network`](@ref).

```@example doc
using Random #hide
Random.seed!(123) #hide
using EcoNetPlot

A = [
    0 0 0 0 0
    1 0 0 0 0
    1 0 0 0 0
    0 1 0 0 0 # Predator eating herbivore 2.
    0 0 1 0 0 # Predator eating herbivore 3.
]
plot_network(A)
```

We see that species are positioned vertically according to their trophic levels.
This is a very simple network, used simply to illustrate
the main goal of the package.
But plotting can be performed on more complex networks.
For example, we can add predators feeding on herbivores in the previous network.

```@example doc
```
Interestingly, we can see that two energy channels are well-separated
thanks to the UMAP embedding.

### From a `Foodweb` or a `Model`

You can plot a network directly from a `Foodweb` object
of `EcologicalNetworksDynamics`. For example


```@example doc
using EcologicalNetworksDynamics

fw = Foodweb(:niche; S = 10, C = 0.1, reject_cycles = true)
plot_network(fw)
```


Not that non-trophic interactions are also supported.

```@example doc
fw = Foodweb(:niche; S = 15, C = 0.1, reject_cycles = true)
nti = NontrophicLayers(;
    :facilitation => (; C = 0.1),
    :interference => (; C = 0.1),
    :refuge => (; C = 0.1),
    :competition => (; C = 0.1),
)
m = default_model(fw, nti)
plot_network(m)
```

Here is the correspondence of the different colours:
- green: facilitation for recruitment
- blue: refuge provisioning
- pink: interference between predators
- light red: competition for space.

For more details please see the [documentation](https://econetoolbox.github.io/EcologicalNetworksDynamics.jl/) of EcologicalNetworksDynamics.

### Plotting the chilean web

Let's plot an empirical network.
The chilean web from [Kéfi et al., Ecology (2015)](https://doi.org/10.1890/13-1424.1)
is directly available within our package with

```@example doc
chilean_web = EcoNetPlot.get_chilean_web()
```


Now let's plot a more complex figure.
We will use the same layout given by [`get_layout`](@ref),
and plot each interaction type in a separate subplot.


```@example doc
using CairoMakie

fig = Figure(; size = (1_000, 400));
ax1 = Axis(fig[1, 1], title = "trophic") 
ax2 = Axis(fig[1, 2], title = "non-trophic positive") 
ax3 = Axis(fig[1, 3], title = "non-trophic negative") 
custom_layout = get_layout(chilean_web[:trophic])
plot_network!(ax1, chilean_web[:trophic]; layout = :custom, custom_layout)
plot_network!(ax2, chilean_web[:positive]; layout = :custom, custom_layout, edge_color = :lightgreen)
plot_network!(ax3, chilean_web[:negative]; layout = :custom, custom_layout, edge_color = :salmon)
fig
```
