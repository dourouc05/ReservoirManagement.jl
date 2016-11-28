["""
Defines the data structures to define a hydrological basin, with its reservoirs and rivers.
"""]

# Basic data structures, upon which the main river/reservoir/basin data structures are based.
include("bathymetry.jl")
include("hydropower.jl")
include("purpose.jl")
include("output.jl")

# Drainage basin
include("river.jl")
include("river_natural.jl")
include("river_diverted.jl")
include("reservoir.jl")
include("basin.jl")

# Nicer functions for the user, but crosses the boundaries of the previous modules, so must be loaded
# after the other parts.
include("bathymetry_fwd.jl")
include("output_fwd.jl")

# Miscellaneous functions.
include("misc.jl")
