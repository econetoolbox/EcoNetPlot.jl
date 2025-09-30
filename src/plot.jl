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
function plot_network!(
    ax,
    A::AbstractMatrix;
    layout = :umap,
    custom_layout = nothing,
    tl_axis = :y,
    A_facilitation = nothing,
    A_competition = nothing,
    A_refuge = nothing,
    A_interference = nothing,
    edge_width = 1,
    node_size = 30 / 2log10(size(A, 1)),
    node_color = :slategray,
    arrow_size = 10,
    edge_color = :sienna,
    curve_distance = 0.01,
    kwargs...,
)
    S = size(A, 1)
    S > 20 && (arrow_size = 0)
    if layout == :custom
        points = custom_layout
    else
        points = get_layout(A; layout, tl_axis)
    end
    g = SimpleDiGraph(A)
    ilabels = S > 20 ? ["" for i in 1:S] : ["$i" for i in 1:S]
    A_list = [A, A_facilitation, A_competition, A_refuge, A_interference]
    A_list = [x for x in A_list if !isnothing(x)]
    degree = get_degree(A_list...)
    node_size = node_size .* (0.4 .+ 0.8log10.(1 .+ degree))
    graphplot!(
        ax,
        g;
        layout = points,
        edge_width = edge_width / log10(sum(A)),
        node_size,
        arrow_size,
        ilabels,
        edge_color,
        node_color,
        curve_distance_usage = true,
        curve_distance,
        kwargs...,
    )
    A_nti_list = [A_facilitation, A_competition, A_refuge, A_interference]
    color_nti = [:lightgreen, :salmon, :lightblue, :pink]
    for (A_nti, color) in zip(A_nti_list, color_nti)
        if !isnothing(A_nti)
            g_nti = DiGraph(A_nti)
            graphplot!(
                ax,
                g_nti;
                edge_color = color,
                node_color,
                layout = points,
                node_size,
                edge_width = edge_width / log10(sum(A_nti)),
                arrow_size,
                ilabels,
                curve_distance_usage = true,
                curve_distance = 2curve_distance,
            )
        end
    end
    hidexdecorations!(ax)
    hidespines!(ax)
    # fig
end
export plot_network!

function plot_network(A; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_network!(ax, A; kwargs...)
    fig
end

function get_layout(A; layout = :umap, tl_axis = :y)
    S = size(A, 1)
    tl = trophic_levels(A)
    A[diagind(A)] .= 0
    if layout == :random
        embedding = rand(S)
    elseif layout == :umap
        min_dist = 10.0 # Control spacing between points.
        distance = get_path_distance(A)
        n_neighbors = min(3, S - 1)
        embedding = umap(distance, 2; metric = :precomputed, n_neighbors, min_dist)[1, :]
        embedding .-= mean(embedding)
        embedding ./= var(embedding)
    else
        @error "layout should be either :umap or :random."
    end
    points = [Point(x, y) for (x, y) in zip(embedding, tl)]
    points = spread_points(points) # Move overlapping points.
    if tl_axis == :x
        points = [Point(p[2], p[1]) for p in points]
    elseif tl_axis != :y
        @error "tl_axis should be either :x or :y."
    end
    points
end
export get_layout

function get_degree(args...)
    A_tot = sum(args)
    degree = vec(sum(A_tot; dims = 1)) .+ vec(sum(A_tot; dims = 2)) / sum(A_tot)
    degree
end
export get_degree

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
function plot_network(m::END.Model; kwargs...)
    A_dict = Dict()
    nti_list = setdiff(m.topology.edge_types_labels, [:trophic])
    for nti in nti_list
        A_dict[nti] = get_matrix(m, nti)
    end
    A_facilitation = haskey(A_dict, :facilitation) ? A_dict[:facilitation] : nothing
    A_interference = haskey(A_dict, :interference) ? A_dict[:interference] : nothing
    A_refuge = haskey(A_dict, :refuge) ? A_dict[:refuge] : nothing
    A_competition = haskey(A_dict, :competition) ? A_dict[:competition] : nothing
    plot_network(
        Matrix(m.A);
        A_facilitation,
        A_interference,
        A_refuge,
        A_competition,
        kwargs...,
    )
end
export plot_network

function get_matrix(m, label)
    if label == :producers_competition
        res = m.producers.competition.matrix |> collect
    else
        res = getproperty(m, label).links.matrix
    end
    res
end

"""
    spread_points(vec; min_dist = 0.015)

Space points so they do not overlap.
"""
function spread_points(vec; min_dist = 0.015)
    x = [p[1] for p in vec]
    tl = [p[2] for p in vec]
    tl_min, tl_max = extrema(tl)
    tl_range = collect(tl_min:3min_dist:tl_max)
    n = length(tl_range)
    for i in 1:n-1
        in_interval(x) = tl_range[i] <= x < tl_range[i+1]
        tl_ind = findall(in_interval, tl)
        isempty(tl_ind) && continue
        D = distance_matrix(x[tl_ind])
        while minimum(D) < min_dist
            idx_list = findall(<(min_dist), D)
            for idx in idx_list
                i, j = Tuple(idx)
                if x[tl_ind[i]] <= x[tl_ind[j]]
                    x[tl_ind[i]] -= min_dist / 4
                    x[tl_ind[j]] += min_dist / 4
                end
            end
            D = distance_matrix(x[tl_ind])
        end
    end
    [Point(x_val, tl_val) for (x_val, tl_val) in zip(x, tl)]
end

function distance_matrix(x)
    n = length(x)
    D = zeros(n, n)
    for i in 1:n, j in 1:n
        if i == j
            D[i, j] = Inf
        else
            D[i, j] = abs(x[i] - x[j])
        end
    end
    D
end
distance_matrix([1, 2, 3])


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

