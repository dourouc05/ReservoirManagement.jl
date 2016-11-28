["""
From given data, handle statistical operations: computing confidence intervals, and generating new
scenarios from given data.
"""]



include("utils.jl")

include("adam2015.jl")
include("tstudent.jl")



## Confidence intervals.

"""
Compute confidence intervals for the input realisations of the time series.

  * `matrix`: one scenario per year. Hypothesis: scenarios are only *one*-year-long.
  * `level`: the confidence level at which the confidence interval should be computed.
  * `ci_period`: on which the confidence intervals must be computed. Only recognised periods:
    * `week` (52 elements per line in `matrix`),
    * `day` (365 elements per line in `matrix`).
  * `method`: the method to use to compute the confidence intervals:
    * `:TStudent` (default): a standard statistical method is used (t-Student hypothesis testing).
      Those confidence intervals are computed for each step independently, with no correlation considered
      between the time steps,
    * `Adam2015`: a physically-based statistical method, based on a sequence of probability distributions.
      This implementation is based on http://hdl.handle.net/2268/142572 without the flooding part.

It is impossible to compute intervals on time scales smaller than the input data (i.e., if `ci_period` is `day`,
but the data is only available weekly). The output will always be at the same discretisation step as the input
scenarios, which may involve copying of the confidence interval values.

TODO:

  * implement other methods, including with correlation.
  * for Adam2015: implement flooding. (Should impact the mean and high part of the CI, not the low?)
  * replace the hard-coded part with a hash table. Keep hard-coded doc for built-in generators, though they are added by "standard" ways (such as addgenerator!() function).
"""
function computeCI_sub(r::River, level::Float64=.95; ci_period::Period=Week(1), method::Symbol=:TStudent)
  TS = r.scenarios
  period = getPeriod(r)
  if period == Week(1) && ci_period == Day(1)
    error("Cannot compute confidence intervals on inexistent data.")
  end

  if period != Week(1) && period != Day(1)
    error("Unrecognised time series period. Only allowable: `week` and `day`.")
  end
  if ci_period != Week(1) && ci_period != Day(1)
    error("Unrecognised confidence interval period. Only allowable: `week` and `day`.")
  end

  # Merge time series as required and make a matrix of it.
  if period != ci_period
    period_f = (period == Day(1)) ? day : (period == Week(1)) ? week : error("Unrecognised period.")
    TS = [collapse(ts, values -> sum(values), period=period_f) for ts in TS] # TODO: test this line hard. Doesn't work properly if period == ci_period (last value much higher than the others).
  end

  # Actual computation.
  # TODO: replace this with a hash table symbol <-> function to make it easier to extend? (Even outside this library!)
  if method == :TStudent
    return computeCI_TStudent(r, level, ci_period=ci_period)
  elseif method == :Adam2015
    return computeCI_Adam2015(r, level)
  end
end

function computeCI(r::River, level::Float64; kwargs...)
  # Sizes? means: (T,). confidence_intervals: (T, 2).
  means, ci_low, ci_high = computeCI_sub(r, level; kwargs...)

  r_mean = copy(r, scenarios=[means])
  r_low  = copy(r, scenarios=[ci_low])
  r_high = copy(r, scenarios=[ci_high])
  return r_mean, r_low, r_high
end

function computeCI(reservoir::Reservoir, level::Float64; ci_period::Period=Week, method::Symbol=:TStudent)
  lnr = [computeCI(r, level, ci_period=ci_period, method=method) for r in naturalRivers(reservoir)]
  ldr = [computeCI(r, level, ci_period=ci_period, method=method) for r in divertedRivers(reservoir)]

  rivers_mean = NaturalRiver[r[1] for r in lnr]
  rivers_low  = NaturalRiver[r[2] for r in lnr]
  rivers_high = NaturalRiver[r[3] for r in lnr]
  divertedRivers_mean = DivertedRiver[r[1] for r in ldr]
  divertedRivers_low  = DivertedRiver[r[2] for r in ldr]
  divertedRivers_high = DivertedRiver[r[3] for r in ldr]

  variantBase = ""
  if reservoir.variant != "vanilla"
    variantBase = reservoir.variant * " "
  end
  variantBase *= "CI@" * string(level * 100) * "% "

  # TODO: add empirical method? (95% of scenarios within bounds, in other words: no real statistics.)
  if method == :TStudent
    techniqueSuffix = "t-Student"
  elseif method == :Adam2015
    techniqueSuffix = "Adam 2015"
  else
    error("Unknown CI technique.")
  end

  reservoir_mean = copy(reservoir, variant=variantBase * "(mean)" * ", " * techniqueSuffix, rivers_in=rivers_mean, rivers_diverted=divertedRivers_mean)
  reservoir_low  = copy(reservoir, variant=variantBase *  "(low)" * ", " * techniqueSuffix, rivers_in=rivers_low,  rivers_diverted=divertedRivers_low)
  reservoir_high = copy(reservoir, variant=variantBase * "(high)" * ", " * techniqueSuffix, rivers_in=rivers_high, rivers_diverted=divertedRivers_high)
  return reservoir_mean, reservoir_low, reservoir_high
end



## Scenario generation from statistics.
"""
Generate a new `River` object based on the given `river`, with its scenarios replaced by newly generated ones.
The generation is parameterised by a number of scenarios to generate `n`, a periodicity for the output `period`,
a generation `method`, and a scenario duration `scenarioDuration` (which must be in years). The parameter `iterator`
can be set to `true` to have an iterator over scenarios (a full reservoir or river being computed at a time) when
multiple ones are requested.

The following methods are implemented:

  * `Adam2015`: a physically-based statistical method, based on a sequence of probability distributions.



Implement new scenario generation functions
===========================================

Mandatory interface to make a new scenario generation function:
    `f(r::River, period::Period, currentYear::Int; kwargs...)`
This function then generates a new scenario for the given `river` for one hydrological year (`currentYear` is only used
to produce good-looking time series objects). The generated data should have a periodicity of `period`, i.e. a value
is output each `period` from the beginning of the hydrological year. The generator might use more parameters than those
ones: they are given as keyword arguments (`kwargs`) through this function.

Steps:

  * create a function with the previous prototype
  * decide for a new keyword (such as :Adam2015); if the same statistical principles are used elsewhere,
    the same symbol should be used
  * modify the function `generateScenarios` to recognise the new keyword and associate the new function
  * update the documentation above

TODO: replace the hard-coded part with a hash table (store the function and the suffix in the same place). Keep hard-coded doc for built-in generators, though they are added by "standard" ways (such as addgenerator!() function).
"""
immutable RiverScenarioIterator
  r::River
  n::Int
  period::Period
  method::Symbol
  scenarioDuration::Period
  kwargs
end
Base.start(::RiverScenarioIterator) = 1
Base.next(i::RiverScenarioIterator, s::Int) = (generateScenarios(i.r, 1, i.period, i.method, i.scenarioDuration; i.kwargs...), s + 1)
Base.done(i::RiverScenarioIterator, s::Int) = s >= i.n
Base.eltype(i::RiverScenarioIterator) = River
Base.length(i::RiverScenarioIterator) = i.n
eachscenario(r::RiverScenarioIterator) = r

function generateScenarios(r::River, n::Int=1, period::Period=Day(1), method::Symbol=:Adam2015, scenarioDuration::Period=Year(1); kwargs...)
  @assert ! Dates.periodisless(scenarioDuration, Year(1)) "Scenario generation only allows generation for time periods that are multiple of years."

  if method == :Adam2015
    f = generateScenarios_Adam2015
  else
    error("Unrecognised method: " * string(method))
  end

  hasIterator = length(filter(e -> e[1] == :iterator, kwargs)) > 0 && filter(e -> e[1] == :iterator, kwargs)[1][2]
  hasIterator && deleteat!(kwargs, find(e -> e[1] == :iterator, kwargs))
  if n == 1 || hasIterator
    # No iterator here.
    years_per_scenario = Int(floor(Dates.days(scenarioDuration) / Dates.days(Year(1))))
    total_n = n * years_per_scenario
    generated = TimeSeries.TimeArray{Float64,1,DateTime,Array{Float64,1}}[let pair = f(r, period, y; kwargs...)
                    TimeArray(pair[1], pair[2], colnames(r.scenarios[1]), meta(r.scenarios[1]))
                  end for y in collect(NonLeapYearIterator(1, total_n))]
    return scenarioGeneration(copy(r, scenarios=generated), :Concatenate, years_per_scenario)
  else
    # Provide an iterator.
    return RiverScenarioIterator(r, n, period, method, scenarioDuration, kwargs)
  end
end

immutable ReservoirScenarioIterator
  reservoir::Reservoir
  n::Int
  period::Period
  method::Symbol
  scenarioDuration::Period
  kwargs
end
Base.start(::ReservoirScenarioIterator) = 1
Base.next(i::ReservoirScenarioIterator, s::Int) = (generateScenarios(i.reservoir, 1, i.period, i.method, i.scenarioDuration; i.kwargs...), s + 1)
Base.done(i::ReservoirScenarioIterator, s::Int) = s >= i.n
Base.eltype(i::ReservoirScenarioIterator) = Reservoir
Base.length(i::ReservoirScenarioIterator) = i.n
eachscenario(r::ReservoirScenarioIterator) = r

function generateScenarios(reservoir::Reservoir, n::Int=1, period::Period=Day(1), method::Symbol=:Adam2015, scenarioDuration::Period=Year(1); kwargs...)
  variant = ""
  if reservoir.variant != "vanilla"
    variant = reservoir.variant * " "
  end
  variant *= "generated scenarios"

  if method == :Adam2015
    techniqueSuffix = "Adam 2015"
  else
    error("Unknown scenario generation technique.")
  end
  variant *= ", " * techniqueSuffix

  hasIterator = length(filter(e -> e[1] == :iterator, kwargs)) > 0 && filter(e -> e[1] == :iterator, kwargs)[1][2]
  if n == 1 || hasIterator
    # No iterator here.
    lnr = NaturalRiver[generateScenarios(r, n, period, method, scenarioDuration; kwargs...) for r in naturalRivers(reservoir)]
    ldr = DivertedRiver[generateScenarios(r, n, period, method, scenarioDuration; kwargs...) for r in divertedRivers(reservoir)]
    return copy(reservoir, variant=variant, rivers_in=lnr, rivers_diverted=ldr)
  else
    # Provide an iterator.
    return ReservoirScenarioIterator(reservoir, n, period, method, scenarioDuration, kwargs)
  end
end
