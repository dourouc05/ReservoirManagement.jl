module ReservoirManagement

  using TimeSeries
  using Base.Dates
  using DataStructures
  using SIUnits
  using SIUnits.ShortUnits
  
  using JuMP
  using Gurobi
  
  using StatsBase
  using HypothesisTests
  using Distributions
  using ExtremeValueDistributions # Strange dependency (not published in METADATA), try to reimplement somewhere
  # (and make it fully compatible with Julia 0.4, without warnings).
  
  using DataFrames
  using Colors
  using Plots
  using StatPlots
  
  # The generalised distributions are defined in both Distributions and ExtremeValueDistributions, but those are identical.
  GeneralizedExtremeValue = Distributions.GeneralizedExtremeValue
  GeneralizedPareto = Distributions.GeneralizedPareto
  
  # Conflict between standard Dates and SIUnits.
  Second = Dates.Second
  
  # To allow overriding those functions.
  import Base: copy, isempty, start, next, done, eltype, length
  
  # Order of inclusion is important!
  include("utils.jl")
  include("reservoir/_ds.jl")
  include("read.jl")
  include("stats/_stats.jl")
  include("scenariosgeneration.jl")
  include("solvers/_solvers.jl")
  include("plots_dataframe.jl")
  include("plots_plots.jl")
  include("report.jl")

end 
