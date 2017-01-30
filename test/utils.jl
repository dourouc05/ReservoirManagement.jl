facts("Utilities") do
  context("Units") do
    nb = 1.
    h = 1.m
    v = 1.m^3
    vv = (v, v)
    q = 1.m^3/s
    p = 1.W

    @fact h::Height --> h
    @fact v::Volume --> v
    @fact q::Discharge --> q
    @fact p::Power --> p

    @fact ReservoirManagement._from_unitful(h)::Float64 --> ReservoirManagement._from_unitful(h)
    @fact ReservoirManagement._from_unitful(v)::Float64 --> ReservoirManagement._from_unitful(v)
    @fact ReservoirManagement._from_unitful(q)::Float64 --> ReservoirManagement._from_unitful(q)
    @fact ReservoirManagement._from_unitful(p)::Float64 --> ReservoirManagement._from_unitful(p)
    @fact ReservoirManagement._from_unitful(nb)::Float64 --> ReservoirManagement._from_unitful(nb)
    @fact ReservoirManagement._from_unitful(vv)::Tuple{Float64, Float64} --> ReservoirManagement._from_unitful(vv)
  end

  context("Season separation: SeasonsSeparator") do
    f = y -> DateTime(y, October, 1)
    @fact_throws SeasonsSeparator(Dict(:a => f), Dict(:b => f)) # IDs mismatch
    @fact typeof(SeasonsSeparator(Dict(:a => f), Dict(:a => f))) <: SeasonsSeparator --> true # No error.
  end

  context("Season separation: getters") do
    @fact length(seasons(HydrologicalWesternEurope)) --> 2
    @fact in(:Dry, seasons(HydrologicalWesternEurope)) --> true
    @fact in(:Wet, seasons(HydrologicalWesternEurope)) --> true

    @fact isSeasonKnown(HydrologicalWesternEurope, :Strange) --> false
    @fact isSeasonKnown(HydrologicalWesternEurope, :Dry) --> true
    @fact isSeasonKnown(HydrologicalWesternEurope, :Wet) --> true

    @fact_throws ensureSeasonKnown(HydrologicalWesternEurope, :Strange)
    @fact ensureSeasonKnown(HydrologicalWesternEurope, :Dry)::Void --> ensureSeasonKnown(HydrologicalWesternEurope, :Dry)
    @fact ensureSeasonKnown(HydrologicalWesternEurope, :Wet)::Void --> ensureSeasonKnown(HydrologicalWesternEurope, :Wet)
  end

  context("Season separation: HydrologicalWesternEurope and getters") do
    @fact begins(HydrologicalWesternEurope, :Wet, 2000) --> Date(2000, October, 1)
    @fact ends(HydrologicalWesternEurope, :Wet, 2000) --> Date(2001, April, 30)
    @fact begins(HydrologicalWesternEurope, :Dry, 2000) --> Date(2001, May, 1)
    @fact ends(HydrologicalWesternEurope, :Dry, 2000) --> Date(2001, September, 30)

    @fact length(HydrologicalWesternEurope, :Wet, 2000) --> Day(212)
    @fact length(HydrologicalWesternEurope, :Dry, 2000) --> Day(153)
    @fact length(HydrologicalWesternEurope, :Wet, 2000, Week(1)) --> 30
    @fact length(HydrologicalWesternEurope, :Dry, 2000, Week(1)) --> 22

    bs, es = seasonBoundariesGuarantee(HydrologicalWesternEurope, Week(1), 2000)
    @fact length(bs[:Wet] : Week(1) : es[:Wet]) --> 30
    @fact length(bs[:Dry] : Week(1) : es[:Dry]) --> 22
    @fact length(bs[:Wet] : Week(1) : es[:Dry]) --> 52
  end

  context("Season separation: separation for HydrologicalWesternEurope, no special guarantee") do
    ts = gen_ts(25 * 52, 2000)
    new_dates = [DateTime(2000, January, 1) + Week(i) for i in 0:length(ts.timestamp) - 1]
    ts = TimeArray(new_dates, values(ts), colnames(ts), meta(ts))

    ts_sep = separate(ts, HydrologicalWesternEurope, 2000)
    ts1w = ts_sep[:Wet]
    ts1w_ts = timestamp(ts1w)
    ts1d = ts_sep[:Dry]
    ts1d_ts = timestamp(ts1d)

    @fact ts1w --> not(ts1d)
    for t in ts1w_ts; @fact in(t, ts1d_ts) --> false; end
    for t in ts1d_ts; @fact in(t, ts1w_ts) --> false; end
    for t in ts1w_ts; @fact t >= begins(HydrologicalWesternEurope, :Wet, 2000) --> true; end
    for t in ts1w_ts; @fact t <=   ends(HydrologicalWesternEurope, :Wet, 2000) --> true; end
    for t in ts1d_ts; @fact t >= begins(HydrologicalWesternEurope, :Dry, 2000) --> true; end
    for t in ts1d_ts; @fact t <=   ends(HydrologicalWesternEurope, :Dry, 2000) --> true; end

    # This separator does not have the property that all years are 52-week long.
    for y in 2000:2024
      ts_sep = separate(ts, HydrologicalWesternEurope, 2000)
      @fact length(ts_sep[:Wet]) + length(ts_sep[:Dry]) --> less_than_or_equal(53)
    end
  end

  context("Season separation: separation for HydrologicalWesternEurope, with 52-week year guarantee") do
    ts = gen_ts(25 * 52, 2000)
    new_dates = [DateTime(2000, January, 1) + Week(i) for i in 0:length(ts.timestamp) - 1]
    ts = TimeArray(new_dates, values(ts), colnames(ts), meta(ts))

    ts_sep = separate(ts, HydrologicalWesternEurope, 2000, guaranteeTimeStep=true)
    ts1w = ts_sep[:Wet]
    ts1w_ts = timestamp(ts1w)
    ts1d = ts_sep[:Dry]
    ts1d_ts = timestamp(ts1d)

    @fact ts1w --> not(ts1d)
    for t in ts1w_ts; @fact in(t, ts1d_ts) --> false; end
    for t in ts1d_ts; @fact in(t, ts1w_ts) --> false; end
    for t in ts1w_ts; @fact t >= begins(HydrologicalWesternEurope, :Wet, 2000) - Week(1) --> true; end
    for t in ts1w_ts; @fact t <=   ends(HydrologicalWesternEurope, :Wet, 2000) + Week(1) --> true; end
    for t in ts1d_ts; @fact t >= begins(HydrologicalWesternEurope, :Dry, 2000) - Week(1) --> true; end
    for t in ts1d_ts; @fact t <=   ends(HydrologicalWesternEurope, :Dry, 2000) + Week(1) --> true; end

    # Specific property of HydrologicalWesternEurope: each year is 52-week long, exactly.
    for y in 2000:2024
      ts_sep = separate(ts, HydrologicalWesternEurope, 2000)
      @fact length(ts_sep[:Wet]) + length(ts_sep[:Dry]) --> 52
    end
  end

  context("Non leap years iterator") do
    for y in collect(NonLeapYearIterator(1900, 500)) # Strange rules every 400 years, hence go for 500.
      @fact y --> not(isleapyear)
    end
  end

  context("Making scenarios from read time series") do
    # Generate some hourly data for a few years.
    hours_per_year = [24 * daysinyear(y) for y in 2000:2005]
    vals = [2 + sin(2π * t / 52) + rand() for t in 1:sum(hours_per_year)]
    dats = [DateTime(2000, January, 1) + Hour(i) for i in 0:sum(hours_per_year) - 1]
    TS = TimeArray(dats, vals, ["Discharge (m³/s)"], Hour(1)) # m^3/s, recorded every hour for six calendar years.

    # Check these randomly generated inputs are consistent.
    @fact timestamp(TS) --> unique(timestamp(TS))
    @fact sum(values(TS) .< 0.) --> 0

    r = makeScenarios(TS, scenarioDuration=Year(1), sampleDuration=Day(1), resamplingRate=1,
                          yearBeginning=begins(HydrologicalWesternEurope, :Wet, 2000))

    # Check the results for r. Complex logic for the values, so quite limited for now...
    @fact length(r) --> 5

    first_date = begins(HydrologicalWesternEurope, :Wet, 2000)
    first_value_origs = TS[first_date : Hour(1) : first_date + Hour(23) + Minute(59)]
    @fact values(r[1])[1] --> roughly(sum(values(first_value_origs) * 3600 / 1.e6))
    for i in 1:5
      y = 2000 + i - 1
      beg_date = begins(HydrologicalWesternEurope, :Wet, y)
      end_date = ends(HydrologicalWesternEurope, :Dry, y)
      if isleapyear(y + 1) # The leap day is in the calendar year after y (always in February,
        # while the hydrological year starts in October).
        end_date -= Day(1)
      end

      @fact meta(r[i]) --> Day(1)
      @fact timestamp(r[i])[1] --> beg_date
      @fact timestamp(r[i])[end] --> end_date
    end
  end

  context("Ensuring an array has the right size, with integer copies allowed to reach that size") do
    array = collect(1.:5.)

    @fact toLength(array, 5) --> array
    @fact toLength(array, 4) --> array[1:4]
    @fact_throws toLength(array, 6) # To reach 6, would need to copy parts of the array.
    @fact toLength(array, 10) --> vcat(array, array) # Exactly two copies of the array.
  end
end
