project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))

"""
    get_chilean_web()

Return the chilean web from Kéfi et al., Ecology (2015).
Contains three type of interactions: trophic, postive and negative.
The object returned is a dictionary, the adjacency of each interaction type can be accessed with its name.

```julia
chilean_web = EcoNetPlot.get_chilean_web()
chilean_web[:trophic] # Access trophic interactions.
chilean_web[:positive] # Access positive interactions.
chilean_web[:negative] # Access negative interactions.
```

Article DOI: https://doi.org/10.1890/13-1424.1
Data DOI: https://doi.org/10.5061/dryad.b4vg0
"""
function get_chilean_web()
    Dict(
        :trophic => deserialize(project_path("data/A_trophic.bin")),
        :positive => deserialize(project_path("data/A_positive.bin")),
        :negative => deserialize(project_path("data/A_negative.bin")),
    )
end
