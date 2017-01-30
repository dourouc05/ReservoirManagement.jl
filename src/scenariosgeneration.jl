["""
Defines generic functions to generate longer scenarios from a list of existing scenarios.
"""]



function _test_unique_dates_inputs(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1})
  try
    vcat(TS...)
  catch e
    error("Dates are not unique in input scenarios! Some dates are present in more than one scenario.")
  end
end

"""
Copy the input scenarios the given number of times along time. More formally, each scenario `s` is copied `number`
times.

For example, if `number = 2`, then `s` becomes `s s`; if `number = 3`, `s` becomes `s s s`.
"""
function scenarioDuplication(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, number::Int)
  _test_unique_dates_inputs(TS)
  n_scenarios = length(TS)

  # To generate scenario i: take input scenario at index i, and take it `number` of time.
  # This technique needs some sophistication: update the years inside each time series so that all of them are unique.
  return TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}[
            vcat([
                map(
                    (timestamp, values) -> (timestamp + Year(1000 * (j - 1)), values), # Add a few millenia to generate new dates.
                    TS[i]
                ) for j in 1:number
            ]...) for i in 1:n_scenarios
         ]
end

"""
Concatenate `number` scenarios to make longer ones, each output scenario containing `number` consecutive input
scenarios. No scenario is ever used twice.

For example, the six scenarios `s1` to `s6` become the three scenarios `s1 s2`, `s3 s4`, and `s5 s6` if `number = 2`.
"""
function scenarioConcatenation(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, number::Int)
  _test_unique_dates_inputs(TS)
  n_output_scenarios = try Int(length(TS) // number)
  catch e
    error("Given number of copies does not divide the number of available scenarios, cannot concatenate properly: " *
          "with " * string(length(TS)) * " available scenarios, " *
          "to output a total of " * string(number) * " scenarios, " *
          "each output scenario would use approximately " * string(length(TS) / number) * " input scenarios.")
  end

  # To generate scenario i: take input scenario at index i, and take it `number` of time.
  # This technique needs some sophistication: update the years inside each time series so that all of them are unique.
  return TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}[
            vcat([TS[(i - 1) * number + j] for j in 1:number]...) for i in 1:n_output_scenarios
         ]
end

"""
Merge the input scenarios: each output scenario is one scenario then the `number` following ones in the input matrix.
Intra- and inter-annual correlations are kept if input scenarios are ordered with time.

For example, the three scenarios `s1` to `s3` become the three scenarios `s1 s2` and `s2 s3` if `number = 2`.
"""
function scenarioMerging(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, number::Int)
  _test_unique_dates_inputs(TS)
  n_scenarios = length(TS)
  new_n_scenarios = n_scenarios - (number - 1) # Number of scenarios to generate.

  # To generate scenario i: take input scenarios starting at index i, and take `number` of them.
  return TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}[
            vcat(TS[i : i + (number - 1)]...) for i in 1:new_n_scenarios
          ]
end

"""
Mix the input scenarios: each output scenario is the concatenation of `number` scenarios, taken arbitrarily
in the set of scenarios (the same one can be used multiple times to generate a given output scenario).
All such combinations are produced. Intra-annual correlations are kept if input scenarios are ordered
with time. Inter-annual correlations are lost.

For example, the three scenarios `s1` to `s3` become nine scenarios if `number = 2`: `s1 s1`, `s1 s2`, `s1 s3`,
`s2 s1`, `s2 s2`, `s2 s3`, `s3 s1`, `s3 s2`, and `s3 s3`.
"""
function scenarioMixing(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, number::Int)
  _test_unique_dates_inputs(TS)

  # Prepend the elements of elements n times to each other. For list(1, 2) and n=2, generate list(list(1, 1), list(1, 2), list(2, 1), list(2, 2)).
  function prependAll(acc::DataStructures.Cons{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}},
                      elements::DataStructures.Cons{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}},
                      n::Int)
    if n <= 1
      return acc
    end

    # For each element of `acc`, generate a list where each element of `elements` is appended; make a list out of these.
    # Add a few millenia to generate new dates, based on the number of already added time series (number - n), with number the total number of series to add, and n the remaining number.
    new_acc = cat( # Concatenate lists.
                  [ # Each list corresponds to one element of the accumulator, with the elements appended.
                      map(e -> vcat(a, map((timestamp, values) -> (timestamp + Year((number - n + 1) * 1000), values), e)),
                          elements) for a in acc]...)

    # Then go for recursivity.
    return prependAll(new_acc, elements, n - 1)
  end

  # Actual mixing, make an array out of the list (prependAll).
  TS_list = cat([list(TS[i]) for i in 1:length(TS)]...)
  return TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}[sc for sc = prependAll(TS_list, TS_list, number)]
end



scenarioGenerationAlgorithms = Dict{Symbol, Tuple{Function, AbstractString}}()

"""
Registers a scenario generation function. The three arguments are:

  * the symbol that the user will use to call the algorithm (like `:Duplicate` or `:Mix`)
  * the function that implements the scenario generation algorithm. It must take two arguments:

    * an array of time series (each one has the same duration)
    * an integer (`number`)

    It returns an array of time series (whose length is `number` times larger than in the input array) with the algorithm applied.

  * the name of the algorithm
"""
function addScenarioGenerationAlgorithm!(s::Symbol, f::Function, d::AbstractString)
  scenarioGenerationAlgorithms[s] = (f, d)
end

_scgen_symbol_from_function(f::Function) = collect(keys(filter((k, v) -> v[1] == f, scenarioGenerationAlgorithms)))[1]
_scgen_function_from_symbol(s::Symbol) = scenarioGenerationAlgorithms[s][1]
_scgen_string_from_symbol(s::Symbol) = scenarioGenerationAlgorithms[s][2]
_scgen_string_from_function(f::Function) = _scgen_string_from_symbol(_scgen_symbol_from_function(f))
_scgen_string_from_any(a::Symbol) = _scgen_string_from_symbol(a)
_scgen_string_from_any(a::Function) = _scgen_string_from_function(a)

addScenarioGenerationAlgorithm!(:Duplicate,   scenarioDuplication,   "duplicate")
addScenarioGenerationAlgorithm!(:Concatenate, scenarioConcatenation, "concatenate")
addScenarioGenerationAlgorithm!(:Merge,       scenarioMerging,       "merge")
addScenarioGenerationAlgorithm!(:Mix,         scenarioMixing,        "mix")



scenarioGeneration(r::River, f::Function, number::Int) = copy(r, scenarios=f(r.scenarios, number)) # TODO: Is this useful, shouldn't all user code directly use symbols (no more functions)? Looks like half-baked refactoring.

"""
Applies a scenario generation algorithm on a river (natural or diverted).
"""
scenarioGeneration(r::River, s::Symbol, number::Int) = scenarioGeneration(r, _scgen_function_from_symbol(s), number)

"""
Applies a scenario generation algorithm on a reservoir (i.e. all its rivers).
"""
function scenarioGeneration(reservoir::Reservoir, f::Union{Symbol, Function}, number::Int) # TODO: Remove the Function parameter?
  newRivers = NaturalRiver[scenarioGeneration(r, f, number) for r in reservoir.rivers_in]
  newDivertedRivers = DivertedRiver[scenarioGeneration(r, f, number) for r in reservoir.rivers_diverted]

  # Update the variant name only if "vanilla", i.e. default reservoir.
  newVariant = reservoir.variant
  if reservoir.variant == "vanilla"
    newVariant = _scgen_string_from_any(f)
  end

  return copy(reservoir, variant=newVariant, rivers_in=newRivers, rivers_diverted=newDivertedRivers)
end

# TODO: autogenerate this from `addScenarioGenerationAlgorithm!`?
scenarioDuplication(r::Reservoir, number::Int) = scenarioGeneration(r, :Duplicate, number)
scenarioDuplication(r::River, number::Int) = scenarioGeneration(r, :Duplicate, number)
scenarioConcatenation(r::Reservoir, number::Int) = scenarioGeneration(r, :Concatenate, number)
scenarioConcatenation(r::River, number::Int) = scenarioGeneration(r, :Concatenate, number)
scenarioMerging(r::Reservoir, number::Int) = scenarioGeneration(r, :Merge, number)
scenarioMerging(r::River, number::Int) = scenarioGeneration(r, :Merge, number)
scenarioMixing(r::Reservoir, number::Int) = scenarioGeneration(r, :Mix, number)
scenarioMixing(r::River, number::Int) = scenarioGeneration(r, :Mix, number)




"""
Keeps a series of indices for each of the given time series:

  * `firstKept`: keep elements starting at this index
  * `length`: keep this number of elements
  * `lastKept`: keep elements up to this index

Exactly two of these three parameters must be provided.
"""
function scenarioShift(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}; firstKept::Int=-1, length::Int=-1, lastKept::Int=-1)
  if firstKept * length * lastKept >= 0
    error("At most two parameters can be given amont firstKept, length, lastKept")
  end

  if firstKept > 0 && length > 0
    return [ts[firstKept : (firstKept + length)] for ts in TS]
  elseif firstKept > 0 && lastKept > 0
    return [ts[firstKept : lastKept] for ts in TS]
  elseif length > 0 && lastKept > 0
    return [ts[(lastKept - length) : lastKept] for ts in TS]
  end
end

"""
Keeps a series of indices for each of the scenarios of the given river:

  * `firstKept`: keep elements starting at this index
  * `length`: keep this number of elements
  * `lastKept`: keep elements up to this index

Exactly two of these three parameters must be provided.
"""
scenarioShift(r::River; firstKept::Int=-1, length::Int=-1, lastKept::Int=-1) = copy(r, scenarios=scenarioShift(r.scenarios, firstKept=firstKept, length=length, lastKept=lastKept))

"""
Keeps a series of indices for each of the scenarios of each river for the given reservoir:

  * `firstKept`: keep elements starting at this index
  * `length`: keep this number of elements
  * `lastKept`: keep elements up to this index

Exactly two of these three parameters must be provided.
"""
function scenarioShift(reservoir::Reservoir; firstKept::Int=-1, length::Int=-1, lastKept::Int=-1)
  newRivers = NaturalRiver[scenarioShift(r, firstKept=firstKept, length=length, lastKept=lastKept) for r in naturalRivers(reservoir)]
  newDivertedRivers = DivertedRiver[scenarioShift(r, firstKept=firstKept, length=length, lastKept=lastKept) for r in divertedRivers(reservoir)]

  return copy(reservoir, rivers_in=newRivers, rivers_diverted=newDivertedRivers)
end
