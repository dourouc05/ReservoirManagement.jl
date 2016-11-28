"""
Describes a river that is diverted to a reservoir. It is akin to a natural
river, with the exception that the actual flow from this source can be decided,
as long as it respects some physical constraints:

  * at least some environmental flow remaining in the river.
  * at most the flow of the river is diverted.
  * at most the capacity of the pipe is diverted.

With respect to the abstract `River`, it defines two new variables to deal with
these new constraints: `environmental_flow` and `maximum_flow`.
"""
immutable DivertedRiver <: River
  name::AbstractString
  scenarios::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1} # 10^6 m^3 / period, indexed with scenarios.
  ext::Dict{AbstractString, Any}

  environmental_flow::Float64 # 10^6 m^3 / s
  maximum_flow::Float64 # 10^6 m^3 / s

  # Private constructor, without units nor defaut values.
  function DivertedRiver(name::AbstractString,
                scenarios::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1},
                ext::Dict{AbstractString, Any},
                environmental_flow::Float64, maximum_flow::Float64)
    o = new(name, scenarios, ext, environmental_flow, maximum_flow)
    _checkConsistency(o)
    return o
  end
end

# Public constructors, with units and keyword arguments.
DivertedRiver(;name::AbstractString=error("Missing reservoir name"),
               scenarios::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}=error("Missing scenarios"),
               ext::Dict{AbstractString, Any}=Dict{AbstractString, Any}(),
               environmental_flow::Discharge=0m^3/s, maximum_flow::Discharge=(Inf)m^3/s) =
  DivertedRiver(name, scenarios, ext, _from_unitful(environmental_flow), _from_unitful(maximum_flow))

copy(r::DivertedRiver; name::AbstractString=r.name,
                       scenarios::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}=r.scenarios,
                       ext::Dict{AbstractString, Any}=r.ext,
                       environmental_flow::Union{Discharge, Float64}=r.environmental_flow,
                       maximum_flow::Union{Discharge, Float64}=r.maximum_flow
    ) =
  # Always call the internal constructor!
  DivertedRiver(name, scenarios, ext, _from_unitful(environmental_flow), _from_unitful(maximum_flow))

Base.isempty(r::DivertedRiver) = false




getEnvironmentalFlow(r::DivertedRiver) = r.environmental_flow
getMaximumFlow(r::DivertedRiver) = r.maximum_flow
getEnvironmentalFlow(r::DivertedRiver, p::Period) = r.environmental_flow * period2seconds(p) # m^3/s -> m^3/period
getMaximumFlow(r::DivertedRiver, p::Period) = r.maximum_flow * period2seconds(p) # m^3/s -> m^3/period

function getMaximumAllowableFlow(r::DivertedRiver, scenario::Int, t::Int)
  period = getPeriod(r)

  actual_flow = getScenario(r, scenario, t)
  remaining_flow = getEnvironmentalFlow(r, period)
  maximum_diverted = getMaximumFlow(r, period)

  available_flow = actual_flow - remaining_flow
  available_flow = (available_flow > 0.) ? available_flow : 0.
  return (available_flow > maximum_diverted) ? maximum_diverted : available_flow
end

environmentalFlow(r::DivertedRiver) = getEnvironmentalFlow(r)
maxFlow(r::DivertedRiver) = getMaximumFlow(r)
environmentalFlow(r::DivertedRiver, p::Period) = getEnvironmentalFlow(r, p)
maxFlow(r::DivertedRiver, p::Period) = getMaximumFlow(r, p)
maxAllowableFlow(r::DivertedRiver, scenario::Int, t::Int) = getMaximumAllowableFlow(r, scenario, t)
