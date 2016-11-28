"""
Describes a river that naturally flows to a reservoir. Its fields correspond
exactly to those of the abstract `River`.
"""
immutable NaturalRiver <: River
  name::AbstractString
  scenarios::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1} # 10^6 m^3 / period, indexed with scenarios.
  ext::Dict{AbstractString, Any}

    # Private constructor, without units nor defaut values.
  function NaturalRiver(name::AbstractString,
                        scenarios::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1},
                        ext::Dict{AbstractString, Any})
    o = new(name, scenarios, ext)
    _checkConsistency(o)
    return o
  end
end

# Public constructors, with units and keyword arguments.
NaturalRiver(;name::AbstractString=error("Missing reservoir name"),
              scenarios::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}=error("Missing scenarios"),
              ext::Dict{AbstractString, Any}=Dict{AbstractString, Any}()) =
  # Always call the internal constructor!
  NaturalRiver(name, scenarios, ext)

copy(r::NaturalRiver; name::AbstractString=r.name, scenarios::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}=r.scenarios, ext::Dict{AbstractString, Any}=r.ext) =
  # Always call the internal constructor!
  NaturalRiver(name, scenarios, ext)

Base.isempty(r::NaturalRiver) = false
