# Plotting ecological networks

A plotting package for ecological networks (food webs and more) powered by [`GraphMakie`](https://github.com/MakieOrg/GraphMakie.jl) 
and [`UMAP`](https://github.com/dillondaudert/UMAP.jl).
Direct integration with [`EcologicalNetworksDynamics`](https://github.com/econetoolbox/EcologicalNetworksDynamics.jl).

### Plotting from an adjacency matrix

Let's say that you have an adjacency matrix which represent an ecological networks, say, a food web.
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
    0 0 0 # Producer (1).
    1 0 0 # Hervibore (2).
    1 0 0 # Herbivore (3).
]
plot_network(A)
```

We see that species are positioned vertically according to their trophic levels.
This is a very simple network, used simply to illustrate
the main goal of the package.
But plotting can be performed on more complex networks.
For example, we can add predators feeding on herbivores in the previous network.

```@example doc
A = [
    0 0 0 0 0
    1 0 0 0 0
    1 0 0 0 0
    0 1 0 0 0 # Predator eating herbivore 2.
    0 0 1 0 0 # Predator eating herbivore 3.
]
plot_network(A)
```

We can further add a top predator that feed on the herbivore's predators.

```@example doc
A = [
    0 0 0 0 0 0
    1 0 0 0 0 0
    1 0 0 0 0 0
    0 1 0 0 0 0
    0 0 1 0 0 0
    0 0 0 1 1 0 # Adding top predator.
]
plot_network(A)
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
