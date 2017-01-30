using ReservoirManagement
using FactCheck
using SIUnits
using SIUnits.ShortUnits
using TimeSeries
using JuMP

include("_utils.jl")

include("reservoir_ds.jl")
include("reservoir_river.jl")
include("reservoir_reservoir.jl")
include("utils.jl")
include("scenariosgeneration.jl")
include("stats.jl")
include("solvers.jl")

FactCheck.exitstatus()
