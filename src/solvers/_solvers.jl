# Main data structures for the solvers.
include("ds.jl")
include("ds_accessors.jl")
include("none.jl")

# Modelling helpers.
include("_helpers_inflow.jl")
include("_helpers_output.jl")
include("_helpers_bathymetry.jl")
include("_helpers_purpose.jl")

# Evaluation solvers.
include("evaluate_crossvalidate.jl")
include("evaluate_deviation.jl")
include("evaluate_mc_feasibility.jl")
include("evaluate_purposeshortage.jl")
include("evaluate_hydropower.jl")

# Inverse solvers.
include("inverse_functions.jl")

# Rule curve solvers.
include("rulecurve_min_default.jl")
include("rulecurve_min_safe.jl")
