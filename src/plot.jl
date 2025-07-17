trophic_levels(A) = EcologicalNetworksDynamics.Internals.trophic_levels(A)

"""
    plot_network(A; layout = :umap, tl_axis = :y)

Plot network graph from adjacency matrix.
By default, species are positioned according to their trophic levels
along the y-axis (can be changed with `tl_axis` argument).
Their position along the x-axis is determined
using umap embedding, but other options are possible (see `layout`).

### Argument

  - `A`: the adjacency matrix (filled with zeros and ones, or booleans).

### Keyword arguments

  - `layout`: the type of layout, either 1) :umap which uses umap to compute
    the embedding, 2) :random which sets random values,
    or 3) :aligned which sets all values equal (useful for chains).
  - `tl_axis`: set the axis of the trophic levels, either :x or :y.
  - `kwargs`: keyword argument given to the graphplot function of GraphMakie.
"""
function plot_network(A::AbstractMatrix; layout = :umap, tl_axis = :y, kwargs...)
    S = size(A, 1)
    tl = trophic_levels(A)
    if layout == :random
        embedding = rand(S)
    elseif layout == :umap
        distance = get_path_distance(A)
        n_neighbors = min(3, S - 1)
        embedding = umap(distance, 2; metric = :precomputed, n_neighbors)[1, :]
        embedding .-= mean(embedding)
        @info embedding
    else
        @error "layout should be either :umap or :random."
    end
    if tl_axis == :y
        points = [Point(x, y) for (x, y) in zip(embedding, tl)]
    elseif tl_axis == :x
        points = [Point(x, y) for (x, y) in zip(tl, embedding)]
    else
        @error "tl_axis should be either :x or :y."
    end
    tlmin, tlmax = extrema(tl)
    xticks = tlmin:tlmax
    g = SimpleDiGraph(A)
    nlabels = ["$i" for i in 1:S]
    fig, ax, p = graphplot(
        g;
        layout = points,
        edge_width = 1,
        node_size = 10,
        arrow_size = 10,
        nlabels,
        edge_color = :gray,
        kwargs...,
    )
    ax.xticks = xticks
    hidexdecorations!(ax)
    hidespines!(ax)
    fig
end

"""
    plot_network(fw::END.Foodweb; kwargs...)

Plotting directly from a Foodweb object of EcologicalNetworksDynamics.
"""
plot_network(fw::END.Foodweb_.Matrix; kwargs...) =
    plot_network(default_model(fw); kwargs...)

plot_network(fw::END.Foodweb_.Adjacency; kwargs...) =
    plot_network(default_model(fw); kwargs...)

"""
    plot_network(m::END.Model; kwargs...)

Plotting directly from a Model object of EcologicalNetworksDynamics.
"""
plot_network(m::END.Model; kwargs...) = plot_network(m.A; kwargs...)

export plot_network


function find_all_paths(g::DiGraph, max_depth::Int)
    paths = []
    function dfs(path)
        push!(paths, copy(path))
        if length(path) >= max_depth
            return
        end
        for nbr in outneighbors(g, last(path))
            if nbr ∉ path  # avoid cycles
                dfs([path...; nbr])
            end
        end
    end
    for v in vertices(g)
        dfs([v])
    end
    return paths
end

function cooccurrence_matrix(paths, n::Int)
    C = zeros(Float64, n, n)
    for path in paths
        for i in path, j in path
            C[i, j] += 1
        end
    end
    return C
end

function similarity_to_distance(C::Matrix{Float64})
    C = C ./ maximum(C)
    1 .- C
end

function get_path_distance(A)
    n = size(A, 1)
    paths = find_all_paths(DiGraph(A), 4)
    C = cooccurrence_matrix(paths, n)
    similarity_to_distance(C)
end

