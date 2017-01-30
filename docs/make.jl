using Documenter
using ReservoirManagement

makedocs(
    modules = [ReservoirManagement],
    format = :html,
    sitename = "ReservoirManagement.jl",
    authors = "Thibaut Cuvelier.",
    analytics = "UA-67298363-2",
    pages = Any[
        "Home" => "index.md",
        "Scenario generation" => "scenario_generation.md"
    ]
)
