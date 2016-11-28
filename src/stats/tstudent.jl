function computeCI_TStudent(r::River, level::Float64=.95; ci_period::Period=Week(1))
  TS = getScenarios(r)
  period = getPeriod(r)
  time_steps = length(TS[1])
  ci_low = zeros(Float64, time_steps)
  ci_high = zeros(Float64, time_steps)
  means = zeros(Float64, time_steps)
  matrix = vcat([values(ts)' for ts in TS]...) # Indexes: scenario, time step.

  # Compute the confidence interval for these time steps, one at a time.
  for t in 1:size(matrix, 2)
    data = vec(matrix[:, t])
    means[t] = mean(data)
    period_ci = confint(OneSampleTTest(data), 1 - level)
    ci_low[t] = (period_ci[1] < 0) ? 0 : period_ci[1]
    ci_high[t] = (period_ci[2] < 0) ? 0 : period_ci[2]
  end

  # Copy things as required.
  if period != ci_period # Only possible case: get days in, compute on weeks, must output days.
    # Computations on one week, average of seven days.
    means /= 7
    ci_low /= 7
    ci_high /= 7

    means[1:364] = vec(repmat(means[1:52]', 7, 1))
    means[365] = means[364] # Scenarios are one-year long. Ensure a year is 365 days, and not 364 = 7 * 52 days.

    ci_low[1:364] = vec(repmat(ci_low[1:52]', 7))
    ci_high[1:364] = vec(repmat(ci_high[1:52]', 7))
    ci_low[365, :] = ci_low[364, :]
    ci_high[365, :] = ci_high[364, :]
  end

  # Finally, make new time series out of it.
  ts_means   = TimeArray(timestamp(TS[1]), means,   colnames(TS[1]), meta(TS[1]))
  ts_ci_low  = TimeArray(timestamp(TS[1]), ci_low,  colnames(TS[1]), meta(TS[1]))
  ts_ci_high = TimeArray(timestamp(TS[1]), ci_high, colnames(TS[1]), meta(TS[1]))

  return ts_means, ts_ci_low, ts_ci_high
end
