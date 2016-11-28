["""
Defines reading helpers, from a time series object, with one measure per hour expressed in cubic metres per second.

For time series objects: the `meta` field is a `Period` object.
"""]





"""
From a long time series `ts`, produce a series of scenarios. This function works with any kind of input time series,
as it only filters on precise dates. The output scenarios are time series whose values are expressed as

  * `ts`: the input time series, which should contain at least one `sampleDuration`. It *must* contain measurements
    every hour, expressed in m^3/s. No check is performed on this. No imputation is performed if one value is missing.
  * `scenarioDuration`: the duration for a scenario (preferably, in years, to make it useable for the other
    parts of the package)
  * `sampleDuration`: the duration for each sample in an output scenario (preferably, `Day(1)` or `Week(1)`
    to make it useable for the other parts of the package)
  * `resamplingRate`: the number of output days that should be made of input values. For example, with
    `resamplingRate = 1`, every period of `sampleDuration` within `ts` makes one output value;
    with `resamplingRate = 2`, two distinct values are output for every `sampleDuration`, but within distinct
    scenarios.
  * `yearBeginning`: the beginning of each year; by default, this is the beginning of the hydrological year.
    The `year` component of the date is ignored.
"""
function makeScenarios(ts::TimeSeries.TimeArray;
                       scenarioDuration::Period=Year(1), sampleDuration::Period=Day(1), resamplingRate::Int=1,
                       yearBeginning::DateTime=begins(HydrologicalWesternEurope, :Wet, 1))
  if length(ts) < 1
    error("Empty time series given")
  end
  if Second(sampleDuration) / (resamplingRate * Second(Hour(1))) < 1.0
    error("resamplingRate too small for sampleDuration: some output scenarios do not even have one input sample. ")
  end

  ## First gather the input samples in bins of duration sampleDuration.
  # Rework yearBeginning to ignore year, with magic so it all works...
  setYear = (d::DateTime, y::DateTime) -> DateTime(year(y), month(d), day(d), hour(d), minute(d), second(d))
  firstConsidered = setYear(yearBeginning, minimum(timestamp(ts)))

  # Put the elements in a queue (FIFO) so they can be retrieved efficiently, bin per bin, later on.
  # This loop is the slowest one in this function due to the date search (~ 75% of total time; adding elements to the queue is ~ 15%).
  # Tried options:
  #   - Using findfirst/findlast instead of boolean indexing is much slower.
  #   - Only search within [currentIndex:end] is much faster, but requires to stuff the beginning of the boolean array.
  #     Keeping this beginning up to date (false for zero to currentIndex - 1) makes the loop slower
  #     (especially making sure the new values are zero; without this, the algorithm is wrong).
  #   - Directly using a Deque instead of Queue is also (marginally) slower.
  samples = Queue(Tuple{DateTime, Array{Float64, 1}})
  currentIndex = findfirst(x -> x >= firstConsidered, timestamp(ts)) # Only consider samples from the beginning of the year.
  while currentIndex <= length(ts)
    # Select [first; first + sampleDuration) (i.e. the time step first + sampleDuration is not included,
    # as it belongs to the next period).
    first = timestamp(ts)[currentIndex]
    idx = (timestamp(ts) .>= first) & (timestamp(ts) .< first + sampleDuration)

    enqueue!(samples, (first, values(ts)[idx]))
    currentIndex += sum(idx)
  end

  ## Then make scenarios.
  # First compute the transformation factor between [m^3/s] in ts and [10^6 m^3/sampleDuration] in output.
  factor = Int(Second(sampleDuration)) / (1000000)

  # Then take elements to make scenarios.
  scenarios = Queue(typeof(ts))
  nBinsPerScenario = Int(floor(Dates.days(scenarioDuration) / Dates.days(sampleDuration))) # Also length of a scenario
  while length(samples) >= nBinsPerScenario # Each iteration consumes nBinsPerScenario sample bins.
    ldates = Array{DateTime, 1}(nBinsPerScenario)
    lvalues = zeros(Float64, nBinsPerScenario, resamplingRate)

    # Maybe skip a bin if the first date is too far away from the beginning of the new year.
    firstDate = front(samples)[1]
    if abs(firstDate - setYear(yearBeginning, firstDate)) >= sampleDuration
      dequeue!(samples)
    end

    # Prepare values for each scenario derived from a succession of bins.
    for i in 1:nBinsPerScenario
      bin = dequeue!(samples)
      ldates[i] = bin[1]
      vals = bin[2]

      for j in 1:resamplingRate
        # Consider values from vals, with resamplingRate == 2:
        #       vals    [   1,   2,   3,   4,   5,   6   ]
        #       scenario    1    2    1    2    1    2
        # Divide by the number of values to get an average discharge [m^3/s].
        # Then perform an integration over sampleDuration, and convert units [10^6 m^3/sampleDuration].
        # The integration is constant: consider the discharge takes the value of each input sample during
        #       sampleDuration / resamplingRate.
        # TODO: what about other integration functions?
        subset = vals[j:resamplingRate:end]
        lvalues[i, j] += sum(subset) / length(subset) * factor
      end
    end

    # Ensure timestamps are not too long for the output (no need to store hours when working only with days).
    if sampleDuration >= Day(1)
      ldates = map(d -> trunc(d, Day), ldates)
    end

    # Finally prepare the output by assembling those components.
    for j in 1:resamplingRate
      dates = map(d -> d + (j - 1) * Year(100), ldates) # Avoid two scenarios having the same dates when resampling.
      enqueue!(scenarios, TimeSeries.TimeArray(dates, lvalues[:, j], colnames(ts), sampleDuration))
    end
  end

  # Return the right type (an array, not a stack).
  return collect(TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, scenarios)
end
