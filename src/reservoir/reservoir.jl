"""
Describes a single water reservoir, using volumes (m^3) and discharges (m^3/s).
It is linked to a series of input rivers that provied water, and has to
fulfil some requirements with the output (tap or industry water, environmental
flow, with or without uncertainty).

TODO:

  * what about topology? (where inputs and diversions can be other reservoirs!).
    Store the rivers inside the basin, and only use names here? Could also get
    rid of different arrays for rivers. Could also store information about the links
    to improve the modelling (e.g., add the wave propagation along this edge).

  * consider uncertainty about purposes (like drinking water consumption)?
"""
immutable Reservoir
  name::AbstractString
  variant::AbstractString
  purposes::Array{Purpose, 1}
  hydropower::Array{HydropowerUnit, 1}
  outputs::Array{DamOutput, 1}

  # Reservoir bathymetry.
  capacity::Tuple{Float64, Float64} # 10^6 m^3; minimum is enough to ensure drinking water intake, by hypothesis.
  bathymetry::ReservoirBathymetry

  # Topography.
  rivers_in::Array{NaturalRiver, 1}
  rivers_diverted::Array{DivertedRiver, 1}

  # Private constructor, without units nor defaut values.
  function Reservoir(name::AbstractString, variant::AbstractString, purposes::Array{Purpose, 1},
                     hydropower::Array{HydropowerUnit, 1}, outputs::Array{DamOutput, 1},
                     capacity::Tuple{Float64, Float64}, bathymetry::ReservoirBathymetry,
                     rivers_in::Array{NaturalRiver, 1}, rivers_diverted::Array{DivertedRiver, 1}
                    )
    ## Some consistency check for the inputs.
    if length(rivers_in) + length(rivers_diverted) == 0
      error("Must have at least one input river!")
    end

    ## If there is hydropower, complement the purposes and outputs.
    # First check if these hydropower units are notalready in the given arrays!
    purposeHPs = HydropowerUnit[getHydropowerUnit(p) for p in filter(p -> isHydropower(p), purposes)]
    outputHPs = HydropowerUnit[getHydropowerUnit(o) for o in filter(o -> isHydropower(o), outputs)]
    for hp in hydropower
      if ! in(hp, purposeHPs)
        push!(purposes, DeterministicHydropower(hp))
      end
      if ! in(hp, outputHPs)
        push!(outputs, HydropowerDamOutput(hp))
      end
    end

    ## Build the object.
    obj = new(name, variant, purposes, hydropower, outputs,
              capacity, bathymetry,
              rivers_in, rivers_diverted)

    ## Check consistency.
    # Check consistency of the bathymetry.
    if capacity[1] > capacity[2]
      error("Inconsistent capacities: the minimum is higher than the maximum!")
    end

    # Check consistency of scenarios.
    countScenariosWithConsistency(obj)
    countTimeStepsWithConsistency(obj)
    getPeriodWithConsistency(obj)
    getTimePointsWithConsistency(obj)

    # Check consistency of the ids.
    ids = vec(vcat([getId(p) for p in purposes], [getId(hp) for hp in hydropower], [getId(o) for o in outputs]))
    ids_cnt = filter((k,v) -> v > 1, countmap(ids))
    if length(ids_cnt) > 0
      error("IDs are not consistent! Already met: " * string(ids_cnt) * ".")
    end

    ## Done!
    return obj
  end
end

# Public constructor, with units and keyword arguments.
Reservoir(;name::AbstractString=error("Missing reservoir name"),
           variant::AbstractString="vanilla",
           purposes::Array{Purpose, 1}=Purpose[],
           hydropower::Array{HydropowerUnit, 1}=HydropowerUnit[],
           outputs::Array{DamOutput, 1}=DamOutput[],

           capacity::Tuple{Volume, Volume}=error("Missing reservoir capacities"),
           bathymetry::ReservoirBathymetry=NullReservoirBathymetry(),

           rivers_in::Array{NaturalRiver, 1}=NaturalRiver[],
           rivers_diverted::Array{DivertedRiver, 1}=DivertedRiver[]
          ) = Reservoir(name, variant, purposes, hydropower, outputs,
                        _from_unitful(capacity), bathymetry,
                        rivers_in, rivers_diverted)

copy(r::Reservoir; name::AbstractString=r.name,
                   variant::AbstractString=r.variant,
                   purposes::Array{Purpose, 1}=r.purposes,
                   hydropower::Array{HydropowerUnit, 1}=r.hydropower,
                   outputs::Array{DamOutput, 1}=r.outputs,
                   capacity::Tuple{Union{Float64, Volume}, Union{Float64, Volume}}=r.capacity,
                   bathymetry::ReservoirBathymetry=r.bathymetry,
                   rivers_in::Array{NaturalRiver, 1}=r.rivers_in,
                   rivers_diverted::Array{DivertedRiver, 1}=r.rivers_diverted
    ) =
  Reservoir(name, variant, purposes, hydropower, outputs,
            _from_unitful(capacity), bathymetry,
            rivers_in, rivers_diverted)

Base.isempty(r::Reservoir) = false



### Generic getters.

getName(r::Reservoir) = r.name
getVariant(r::Reservoir) = r.variant
getMinimumCapacity(r::Reservoir) = r.capacity[1]
getMaximumCapacity(r::Reservoir) = r.capacity[2]

getHydropowerUnits(r::Reservoir) = r.hydropower
hasHydropower(r::Reservoir) = length(getHydropowerUnits(r)) > 0
getHydropowerIds(r::Reservoir) = [getId(hp) for hp in getHydropowerUnits(r)]
getPenstockHydropower(r::Reservoir) = sum([hp.max_discharge for hp in getHydropowerUnits(r)])

getPurposes(r::Reservoir) = r.purposes
getOutputs(r::Reservoir) = r.outputs
getBathymetry(r::Reservoir) = r.bathymetry

function _check_purpose_subtype(purposeType::DataType)
  if ! (purposeType <: Purpose)
    error("Wrong type argument: must be a subtype of Purpose; got " * string(purposeType) * ".")
  end
end
function hasPurpose(r::Reservoir, purposeType::DataType)
  _check_purpose_subtype(purposeType)
  return length(filter(t -> isa(t, purposeType), getPurposes(r))) > 0
end
function getPurpose(r::Reservoir, purposeType::DataType)
  _check_purpose_subtype(purposeType)
  # Return an array of purposes, as the type has no unicity guarantee (as opposed to IDs).
  return filter(t -> isa(t, purposeType), getPurposes(r))
end
getPurposeIds(r::Reservoir) = Symbol[getId(o) for o in getPurposes(r)]
hasPurpose(r::Reservoir, s::Symbol) = length(filter(o -> getId(o) == s, getPurposes(r))) > 0
function getPurpose(r::Reservoir, s::Symbol; nullable::Bool=false)
  elt = filter(o -> getId(o) == s, getPurposes(r))
  if nullable
    if length(elt) == 0
      return Nullable{Symbol}()
    else
      return Nullable(elt[1])
    end
  else
    if length(elt) == 0
      error("Purpose " * string(s) * " not found.")
    else
      return elt[1]
    end
  end
end

getOutputIds(r::Reservoir) = Symbol[getId(o) for o in getOutputs(r)]
hasOutput(r::Reservoir, s::Symbol) = length(filter(o -> getId(o) == s, getOutputs(r))) > 0
function getOutput(r::Reservoir, s::Symbol; nullable::Bool=false)
  elt = filter(o -> getId(o) == s, getOutputs(r))
  if nullable
    if length(elt) == 0
      return Nullable{Symbol}()
    else
      return Nullable(elt[1])
    end
  else
    if length(elt) == 0
      error("Purpose " * string(s) * " not found.")
    else
      return elt[1]
    end
  end
end

name(r::Reservoir) = getName(r)
variant(r::Reservoir) = getVariant(r)
minCapacity(r::Reservoir) = getMinimumCapacity(r)
maxCapacity(r::Reservoir) = getMaximumCapacity(r)
penstockHydropower(r::Reservoir) = getPenstockHydropower(r)
hydropowerUnits(r::Reservoir) = getHydropowerUnits(r)
outputs(r::Reservoir) = getOutputs(r)
outputsId(r::Reservoir) = getOutputIds(r)
output(r::Reservoir, s::Symbol; kwargs...) = getOutput(r, s; kwargs...)
bathymetry(r::Reservoir) = getBathymetry(r)
purpose(r::Reservoir, s::Symbol; kwargs...) = getPurpose(r, s; kwargs...)
purpose(r::Reservoir, purposeType::DataType) = getPurpose(r, purposeType)

# Counting things.

"""
Counts the number of scenarios within this reservoir.
"""
countScenarios(r::Reservoir) = countScenarios(first(rivers(r)))

"""
Counts the number of scenarios within this reservoir and ensures the whole data structure is consistent.
"""
function countScenariosWithConsistency(r::Reservoir)
  cnt = -1

  for river in rivers(r)
    if cnt < 0
      cnt = countScenariosWithConsistency(river)
    elseif cnt != countScenariosWithConsistency(river)
      error("Inconsistent data in reservoir " * getName(r) * ": different number of scenarios within the reservoir (while checking river " * getName(river) * ").")
    end
  end

  if cnt < 0
    error("Inconsistent data in reservoir " * getName(r) * ": no time series found in the rivers.")
  end

  return cnt
end

"""
Counts the number of time steps within this reservoir.
"""
countTimeSteps(r::Reservoir) = countTimeSteps(first(rivers(r)))

"""
Counts the number of time steps within this reservoir and ensures the whole data structure is consistent.
"""
function countTimeStepsWithConsistency(r::Reservoir)
  cnt = -1

  for river in rivers(r)
    if cnt < 0
      cnt = countTimeStepsWithConsistency(river)
    elseif cnt != countTimeStepsWithConsistency(river)
      error("Inconsistent data in reservoir " * getName(reservoir) * ": different number of time steps within the reservoir " *
            "(while checking river " * getName(river) * ").")
    end
  end

  if cnt < 0
    error("Inconsistent data in reservoir " * getName(reservoir) * ": no time series found in the rivers.")
  end

  return cnt
end

"""
Retrieves the periodicity of this reservoir data.
"""
getPeriod(r::Reservoir) = getPeriod(getRivers(r)[1])

"""
Retrieves the periodicity of this reservoir data and ensures the whole data structure is consistent.
"""
function getPeriodWithConsistency(r::Reservoir)
  period = :Unknown

  for river in rivers(r)
    if period == :Unknown
      period = getPeriodWithConsistency(river)
    elseif period != getPeriodWithConsistency(river)
      error("Inconsistent data in reservoir " * getName(reservoir) * ": different period within the reservoir " *
            "(while checking river " * getName(river) * ").")
    end
  end

  if period == :Unknown
    error("Inconsistent data in reservoir " * getName(reservoir) * ": no time series found in the rivers.")
  end

  return period::Period
end

"""
Retrieves the time series of this reservoir data, i.e. the time points where historical information is known.
"""
getTimePoints(r::Reservoir) = getTimePoints(getRivers(r)[1])

"""
Retrieves the time instants of this reservoir data, i.e. the time points where historical information is known.
This function also ensures the whole data structure is consistent.
"""
function getTimePointsWithConsistency(r::Reservoir)
  ts = nothing

  for river in rivers(r)
    if ts == nothing
      ts = getTimePointsWithConsistency(river)
    elseif ts != getTimePointsWithConsistency(river)
      error("Inconsistent data in reservoir " * getName(reservoir) * ": different time series within the reservoir " *
            "(while checking river " * getName(river) * ").")
    end
  end

  if ts == nothing
    error("Inconsistent data in reservoir " * getName(r) * ": no time series found in the rivers.")
  end

  return ts::Array{DateTime, 1}
end



### Getters for rivers.

# Iterate through all rivers.
rivers(r::Reservoir) = RiverIterator(r)
naturalRivers(r::Reservoir) = NaturalRiverIterator(r)
divertedRivers(r::Reservoir) = DivertedRiverIterator(r)

getNaturalRivers(r::Reservoir) = r.rivers_in
getDivertedRivers(r::Reservoir) = r.rivers_diverted
getRivers(r::Reservoir) = vcat(getNaturalRivers(r), getDivertedRivers(r))

getNaturalRiver(r::Reservoir, naturalRiver::Int) = r.rivers_in[naturalRiver]
getDivertedRiver(r::Reservoir, divertedRiver::Int) = r.rivers_diverted[divertedRiver]

hasRiver(r::Reservoir, river::AbstractString) = length(filter(nr -> getName(nr) == river, getRivers(r))) >= 1
function getRiver(r::Reservoir, river::AbstractString)
  elt = filter(nr -> getName(nr) == river, getRivers(r))
  if length(elt) < 1
    return Nullable{River}()
  else
    return Nullable(elt[1])
  end
end

# Get river name.
getDivertedRiverNames(reservoir::Reservoir) = UTF8String[r.name for r in divertedRivers(reservoir)]
getNaturalRiverNames(reservoir::Reservoir) = UTF8String[r.name for r in naturalRivers(reservoir)]

function getRiverNames(reservoir::Reservoir)
  names = getNaturalRiverNames(reservoir)
  append!(names, getDivertedRiverNames(reservoir))
  return names
end



### Agregators.

# Total inflow.

"""
Compute the total inflow for this reservoir for the given scenario,
i.e. this function agregates all tributaries to the reservoir.
The output is given along time, either as an array (`asTimeSeries=false`, default),
or as a time series object  (`asTimeSeries=true`).
"""
function totalInflow(r::Reservoir, scenario::Int; asTimeSeries::Bool=false)
  if asTimeSeries
    anyRiver = getRivers(r)[1].scenarios[1]
    return TimeArray(timestamp(anyRiver), totalInflow(r, scenario, asTimeSeries=false), colnames(anyRiver), meta(anyRiver))
  else
    return sum([values(riv.scenarios[scenario]) for riv in eachriver(r)])
  end
end

"""
Compute the total inflow for this reservoir,
i.e. this function agregates all tributaries to the reservoir.
The output is given along time, either as an array (`asTimeSeries=false`, default),
or as a time series object  (`asTimeSeries=true`).
"""
function totalInflow(r::Reservoir; kwargs...)
  if countScenarios(r) > 1
    error("This reservoir has multiple scenarios, please indicate which one to use.")
  end
  return totalInflow(r, 1; kwargs...)
end

"""
Compute the total inflow for this reservoir for the given range of scenarios,
i.e. this function agregates all tributaries to the reservoir.
The output is given along time, either as an array (`asTimeSeries=false`, default),
or as a time series object  (`asTimeSeries=true`).
"""
function totalInflow{T}(r::Reservoir, scenarios::Range{T}; kwargs...)
  kwargs = Dict(kwargs)
  if haskey(kwargs, :asTimeSeries) && kwargs[:asTimeSeries]
    return TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}[totalInflow(r, scenario; kwargs...) for scenario in scenarios]
  else
    return Float64[totalInflow(r, scenario; kwargs...) for scenario in scenarios]
  end
end

# No agregator for outputs or purposes: they have a high variability in the way to take them into account.
# Such a method would only make sense for a subset of all purposes.



### Altering the set of scenarios.

"""
Removes the given set of scenarios from the reservoir, of which a modified copy is given.
"""
function removeScenarios(reservoir::Reservoir, toRemove::IntSet)
  nScenarios = countScenarios(reservoir)
  indices = trues(nScenarios)
  for i in toRemove
    indices[i] = false
  end

  newRivers = NaturalRiver[copy(r, scenarios=r.scenarios[indices]) for r in eachnaturalriver(reservoir)]
  newDivertedRivers = DivertedRiver[copy(r, scenarios=r.scenarios[indices]) for r in eachdivertedriver(reservoir)]

  return copy(reservoir, rivers_in=newRivers, rivers_diverted=newDivertedRivers)
end

removeScenarios(r::Reservoir, toRemove::Int) = removeScenarios(r, IntSet(toRemove))
removeScenarios(r::Reservoir, toRemove::Array{Int, 1}) = removeScenarios(r, IntSet(toRemove))

"""
Limits the scenarios in the reservoir to the given time steps.
"""
function keepTimeSteps(reservoir::Reservoir, timeSteps)
  newRivers = NaturalRiver[keepTimeSteps(r, timeSteps) for r in eachnaturalriver(reservoir)]
  newDivertedRivers = DivertedRiver[keepTimeSteps(r, timeSteps) for r in eachdivertedriver(reservoir)]

  return copy(reservoir, rivers_in=newRivers, rivers_diverted=newDivertedRivers)
end

"""
Limits the scenarios in the river to the given time steps.
"""
keepTimeSteps(r::River, timeSteps) = copy(r, scenarios=[s[timeSteps] for s in r.scenarios])



### River iterators.

# Natural rivers.
immutable NaturalRiverIterator
  r::Reservoir
end

eachnaturalriver(r::Reservoir) = NaturalRiverIterator(r)
Base.start(::NaturalRiverIterator) = 1
Base.next(i::NaturalRiverIterator, s) = (i.r.rivers_in[s], s + 1)
Base.done(i::NaturalRiverIterator, s) = s > length(i.r.rivers_in)
Base.eltype(i::NaturalRiverIterator) = NaturalRiver
Base.length(i::NaturalRiverIterator) = length(i.r.rivers_in)

# Diverted rivers.
immutable DivertedRiverIterator
  r::Reservoir
end

eachdivertedriver(r::Reservoir) = DivertedRiverIterator(r)
Base.start(::DivertedRiverIterator) = 1
Base.next(i::DivertedRiverIterator, s) = (i.r.rivers_diverted[s], s + 1)
Base.done(i::DivertedRiverIterator, s) = s > length(i.r.rivers_diverted)
Base.eltype(i::DivertedRiverIterator) = DivertedRiver
Base.length(i::DivertedRiverIterator) = length(i.r.rivers_diverted)

# Finally for both at once.
# Implementation detail: first natural rivers, then diverted.
# TODO: redo with Iterators.chain()?
immutable RiverIterator
  r::Reservoir
  nIn::Int

  function RiverIterator(reservoir::Reservoir)
    return new(reservoir, length(reservoir.rivers_in))
  end
end

eachriver(r::Reservoir) = RiverIterator(r)
Base.start(::RiverIterator) = 1
Base.next(i::RiverIterator, s) = ((s <= i.nIn) ? i.r.rivers_in[s] : i.r.rivers_diverted[s - i.nIn], s + 1)
Base.done(i::RiverIterator, s) = s > i.nIn + length(i.r.rivers_diverted)
Base.eltype(i::RiverIterator) = River
Base.length(i::RiverIterator) = i.nIn + length(i.r.rivers_diverted)




### Iterate over scenarios in a reservoir.
# Provides a copy of the reservoir with a single scenario.
immutable ScenarioReservoirIterator
  r::Reservoir
  max::Int

  function ScenarioReservoirIterator(reservoir::Reservoir)
    return new(reservoir, countScenarios(reservoir))
  end
end

eachscenario(r::Reservoir) = ScenarioReservoirIterator(r)
Base.start(::ScenarioReservoirIterator) = 1
Base.next(i::ScenarioReservoirIterator, s) = (
  copy(i.r,
       rivers_in=NaturalRiver[copy(riv, scenarios=[riv.scenarios[s]]) for riv in eachnaturalriver(i.r)],
       rivers_diverted=DivertedRiver[copy(riv, scenarios=[riv.scenarios[s]]) for riv in eachdivertedriver(i.r)]),
  s + 1
)
Base.done(i::ScenarioReservoirIterator, s) = s > i.max
Base.eltype(i::ScenarioReservoirIterator) = Reservoir
Base.length(i::ScenarioReservoirIterator) = i.max
