value(x::SIUnits.SIQuantity) = x.val # While not in master: https://github.com/Keno/SIUnits.jl/issues/14

Height = quantity(Float64, SIUnits.Meter)
Volume = quantity(Float64, SIUnits.Meter^3)
Discharge = quantity(Float64, SIUnits.Meter^3/SIUnits.Second)
Power = quantity(Float64, SIUnits.Watt)

_from_unitful(x::Height) = value(x)
_from_unitful(x::Volume) = value(x) / 10^6
_from_unitful(x::Discharge) = value(x) / 10^6
_from_unitful(x::Power) = value(x) / 10^6
_from_unitful(x::Tuple{Volume, Volume}) = (_from_unitful(x[1]), _from_unitful(x[2]))

_from_unitful(x::Float64) = x
_from_unitful(x::Tuple{Float64, Float64}) = x



"""
This set of functions takes as input a yearly time series and output the beginning or the end of the hydrological season.
"""
immutable SeasonsSeparator
  date_begin::Dict{Symbol, Function}
  date_end::Dict{Symbol, Function}

  function SeasonsSeparator(date_begin::Dict{Symbol, Function}, date_end::Dict{Symbol, Function})
    # Check if IDs match. The order does not necessarily match, though! Hence the complex code with looping.
    beg_ids = collect(keys(date_begin))
    end_ids = collect(keys(date_end))

    for id in beg_ids
      if ! in(id, end_ids)
        error("The seasons do not match: found a beginning for season " * string(id) * ", but no end.")
      end
    end
    for id in end_ids
      if ! in(id, beg_ids)
        error("The seasons do not match: found an end for season " * string(id) * ", but no beginning.")
      end
    end

    return new(date_begin, date_end)
  end
end

"""
A seasons separator for Western Europe, based on hydrological years (starts in October, ends in September).
           October   April	|	May			September
             1st     30th   | 1st       30th
                  Wet       |       Dry
"""
HydrologicalWesternEurope = SeasonsSeparator(
  Dict(:Wet => y -> DateTime(y, October, 1),
       :Dry => y -> DateTime(y + 1, May, 1)),
  Dict(:Wet => y -> DateTime(y + 1, April, Dates.daysinmonth(y + 1, April)),
       :Dry => y -> DateTime(y + 1, September, Dates.daysinmonth(y + 1, September)))
)

seasons(seasonsSep::SeasonsSeparator) = collect(keys(seasonsSep.date_begin))
isSeasonKnown(seasonsSep::SeasonsSeparator, season::Symbol) = in(season, seasons(seasonsSep))
function ensureSeasonKnown(seasonsSep::SeasonsSeparator, season::Symbol)
  if ! isSeasonKnown(seasonsSep, season)
    error("Unknown season " * string(seasons) * ". " *
          "The defined seasons for this separator are: " * string(seasons(seasonsSep)))
  end
end
function begins(seasonsSep::SeasonsSeparator, season::Symbol, chosenYear::Int)
  ensureSeasonKnown(seasonsSep, season)
  return seasonsSep.date_begin[season](chosenYear)
end
function ends(seasonsSep::SeasonsSeparator, season::Symbol, chosenYear::Int)
  ensureSeasonKnown(seasonsSep, season)
  return seasonsSep.date_end[season](chosenYear)
end
length(seasonsSep::SeasonsSeparator, season::Symbol, chosenYear::Int) =
  Day(length(begins(seasonsSep, season, chosenYear) : ends(seasonsSep, season, chosenYear)))
function length(seasonsSep::SeasonsSeparator, season::Symbol, chosenYear::Int, period::Period)
  daily = Dates.days(length(seasonsSep, season, chosenYear))
  days_per_period = Dates.days(period)
  return round(Int, daily / days_per_period)
end

"""
Computes the length of each season while guaranteeing that, for the given `period`, all seasons have the same number of
time steps, whatever the value of `chosenYear` is.
"""
function seasonBoundariesGuarantee(seasonsSep::SeasonsSeparator, period::Period, chosenYear::Int)
  ## Idea: start from a basic decomposition based on the actual year; then reassign time steps at the borders so
  ## that each season is closer to its preferred length.

  # Compute the preferred lengths for an arbitrary year (here, 2000, which is not leap).
  preferred = [season => length(seasonsSep, season, 2000, period) for season in seasons(seasonsSep)]

  # Base beginnings and ends, with season lengths.
  bs = [season => begins(seasonsSep, season, chosenYear) for season in seasons(seasonsSep)] # Begin
  es = [season =>   ends(seasonsSep, season, chosenYear) for season in seasons(seasonsSep)] # End
  ls = [season => length(seasonsSep, season, chosenYear, period) for season in seasons(seasonsSep)] # Length

  # Quantify the delta between those values and the preferred ones: if some periods are missing, decide to take them
  # before the beginning of the first season. Otherwise, the delta is not important for now: it will be compensated
  # at the end of the seasons with the next loop.
  total_delta = 0
  for season in seasons(seasonsSep)
    delta = ls[season] - preferred[season]
    if delta < 0 # Length too small for this season: should have started before.
      total_delta += abs(delta)
    elseif delta > 0 # Too long: should have started after, i.e. compensate some shorter seasons.
      total_delta -= abs(delta)
    end
  end
  if total_delta > 0
    season = minimum(bs)[1]
    bs[season] -= total_delta * period
  end

  # Now, adapt the other boundaries in between: start from the first season, then add up the periods up
  # with the reference length.
  ordered_seasons = [pair[2] for pair in sort(collect(zip(values(bs), keys(bs))))]
  es[ordered_seasons[1]] = bs[ordered_seasons[1]] + (preferred[ordered_seasons[1]] - 1) * period
  for i in 2:length(ordered_seasons)
    prev_season = ordered_seasons[i - 1]
    season = ordered_seasons[i]
    bs[season] = es[prev_season] + period
    es[season] = bs[season] + (preferred[season] - 1) * period
  end

  # Done!
  return bs, es
end

"""
Separates seasons for the given time series `ts`. The seasons begin and end as imposed by `seasonsSep`.
The separation is performed for year `chosenYear`. The output is a dictionary, where the time series are indexed
by the season identifier.

When `guaranteeTimeStep` is `true`, then the seasons are such that all years have the same number of
time steps (whose duration is the same as that of `ts`). This option might introduce slight shifts in the season
(i.e. it may begin or end at dates slightly different from `seasonsSep`).
"""
function separate(ts::TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, seasonsSep::SeasonsSeparator,
                  chosenYear::Int=year(minimum(timestamp(ts))); guaranteeTimeStep::Bool=false)
  # Simple case: just do it.
  if ! guaranteeTimeStep
    rets = Dict{Symbol, TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}}()
    for season in seasons(seasonsSep)
      rets[season] = to(from(ts, begins(seasonsSep, season, chosenYear)), ends(seasonsSep, season, chosenYear))
    end
    return rets
  end

  # Make the seasons dictionary with the time step guarantee.
  (bs, es) = seasonBoundariesGuarantee(seasonsSep, meta(ts), chosenYear)
  rets = Dict{Symbol, TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}}()
  for season in seasons(seasonsSep)
    rets[season] = to(from(ts, bs[season]), es[season])
  end
  return rets
end

function separate(ts::TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, seasonsSep::SeasonsSeparator,
                  season::Symbol, chosenYear::Int=year(minimum(timestamp(ts))); guaranteeTimeStep::Bool=false)
  # Simple case: just do it.
  if ! guaranteeTimeStep
    rets = Dict{Symbol, TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}}()
    return to(from(ts, begins(seasonsSep, season, chosenYear)), ends(seasonsSep, season, chosenYear))
  end

  # Make the seasons dictionary with the time step guarantee.
  (bs, es) = seasonBoundariesGuarantee(seasonsSep, meta(ts), chosenYear)
  return to(from(ts, bs[season]), es[season])
end

"""
Separates seasons for the given array of time series.
"""
separate(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, seasons::SeasonsSeparator, season::Symbol, year::Int) =
  [separate(ts, seasons, season, year) for ts in TS]
separate(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, seasons::SeasonsSeparator, season::Symbol) =
  [separate(ts, seasons, season) for ts in TS]
separate(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, seasons::SeasonsSeparator, year::Int) =
  [separate(ts, seasons, year) for ts in TS]
separate(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, seasons::SeasonsSeparator) =
  [separate(ts, seasons) for ts in TS]

"""
Collects seasons for the given time series: throughout the time series, it collects the asked season (e.g., all time steps within the wet season);
for one-year series, it is equivalent to `separate`.
"""
function Base.collect(ts::TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, seasons::SeasonsSeparator, season::Symbol)
  minY = year(minimum(timestamp(ts)))
  maxY = year(maximum(timestamp(ts)))
  return vcat([separate(ts, seasons, season, year) for year in minY:maxY]...)
end



"""
A non-leap year iterator. State: current year in iteration, number of already generated years.

It is more efficient than filtering over years with the built-in function `isleapyear` (Julia 0.4.6):
  @time collect(NonLeapYearIterator(1, 1_000))
       0.000059 seconds (648 allocations: 26.734 KB)
  @time take(filter(x -> ! isleapyear(x), 1:1_500), 1_000)
       0.001695 seconds (1.03 k allocations: 28.041 KB)
"""
immutable NonLeapYearIterator
  start::Int
  n::Int
end

typealias NonLeapYearIteratorState Tuple{Int64, Int64}
Base.start(nly::NonLeapYearIterator) = (nly.start, 0)
function Base.next(nly::NonLeapYearIterator, state::NonLeapYearIteratorState)
  y = state[1] + 1
  if isleapyear(y) # Never two consecutive leap years.
    y = y + 1
  end
  return y, (y, state[2] + 1)
end
Base.done(nly::NonLeapYearIterator, state::NonLeapYearIteratorState) = state[2] >= nly.n



"""
These types allows considering an array of time series as scenarios within the `feasibilityProbability` solver.
"""
immutable ArrayAsIterableReservoirScenarios
  scenarios::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}
end
Base.start(ars::ArrayAsIterableReservoirScenarios) = 1
Base.next(ars::ArrayAsIterableReservoirScenarios, state::Int) =
  ArrayAsIterableReservoirScenario(values(ars.scenarios[state])), state + 1
Base.done(ars::ArrayAsIterableReservoirScenarios, state::Int) = state >= length(ars.scenarios)
eachscenario(lts::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}) =
  ArrayAsIterableReservoirScenarios(lts)

immutable ArrayAsIterableReservoirScenario # Define a new type to avoid overloading totalInflow() with too generic types.
  scenarios::Array{Float64, 1}
end
totalInflow(ars::ArrayAsIterableReservoirScenario) = ars.scenarios



period2seconds(p::Period) = Dates.days(p) * 24 * 60 * 60



"""
Ensures the given `array` has exactly the length `givenLength`, either by removing parts of it, or by copying
the whole array an integer number of times (otherwise, an error is thrown).
"""
function toLength(array::Array{Float64, 1}, givenLength::Int)
  if length(array) == givenLength
    return array
  elseif length(array) > givenLength
    return array[1:givenLength]
  elseif length(array) < givenLength
    n_copies = try
      Int(givenLength // length(array))
    catch
      exact_str = string(givenLength) * "//" * string(length(array)) # Direct string(//) will reduce the fraction.
      error("To match the expected length, would need to copy a noninteger number of times " *
            "(" * exact_str            * " ~= " * string(givenLength / length(array)) * ").")
      #           Exact version (like 1//2)       Approximation (like 0.5)
    end
    return vec(repmat(array, 1, n_copies))
  end
end
