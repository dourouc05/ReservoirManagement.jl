"""
A complete statistical model for one river, with three-stage modelling: yearly average, seasonal averages, daily averages.
"""
immutable Adam2015RiverStatisticalModel
  year::GeneralizedExtremeValue # Annual average; Weibull when xi = -1.
  season_wet::GeneralizedExtremeValue # Difference between seasonal and annual averages, normalised wrt annual average.
  season_dry::GeneralizedExtremeValue
  daily_wet::GeneralizedExtremeValue
  daily_dry::GeneralizedExtremeValue

  seasons_separator::SeasonsSeparator
end

getYearlyAverage(m::Adam2015RiverStatisticalModel) = m.year
function getSeasonAverageDelta(m::Adam2015RiverStatisticalModel, season::Symbol)
  if season == :Wet
    return m.season_wet
  elseif season == :Dry
    return m.season_dry
  else
    error("Unknown season: " * string(season))
  end
end
function getSeasonDailyDelta(m::Adam2015RiverStatisticalModel, season::Symbol)
  if season == :Wet
    return m.daily_wet
  elseif season == :Dry
    return m.daily_dry
  else
    error("Unknown season: " * string(season))
  end
end
getSeasonSeparator(m::Adam2015RiverStatisticalModel) = m.seasons_separator



"""
Based on the data contained in `TS` (one time series per year!), fit a statistical model.

  * `TS`: array of time series (one series per year) to fit the statistical model on.
  * `verbosity`: level of output.
    * `0`: no output.
    * `1`: summary output
    * `2`: detailed output of this algorithm's details
    * `3`: detailed output of this algorithm's and its sub-algorithms' details
    * `4`: including solver output

TODO: implement flooding.
"""
function fit_Adam2015(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1};
                      seasonsSeparator::SeasonsSeparator=HydrologicalWesternEurope,
                      verbosity::Int=1)
  ## Separate every year in two seasons: dry and wet.
  sep_ts = [separate(ts, seasonsSeparator, guaranteeTimeStep=true) for ts in TS]
  wet_TS = [ts[:Wet] for ts in sep_ts]
  dry_TS = [ts[:Dry] for ts in sep_ts]

  # Compute a few useful averages within those seasons.
  ann_avg = [mean(values(ts)) for ts in TS]
  wet_avg = [mean(values(ts)) for ts in wet_TS]
  dry_avg = [mean(values(ts)) for ts in dry_TS]

  if verbosity >= 1
    println(" == ")
    println(" ==== Statistical model fitting on ", length(TS), " years")

    if verbosity >= 2
      println(" ====== Length of wet season: ", length(wet_TS[1]))
      println(" ====== Length of dry season: ", length(dry_TS[1]))
    end
  end

  ## First deal with annual average.
  if verbosity >= 1
    println(" ==== Fitting a generalised extreme values model on the annual average.")
  end

  d_year_params = fit_mle_optim(ExtremeValueDistributions.GeneralizedExtremeValue, ann_avg, [0.5, 0.5, -1], verbose=(verbosity >= 3))
  d_year = Distributions.GeneralizedExtremeValue(d_year_params[1], exp(d_year_params[2]), d_year_params[3]) # ~ <Q>_Y

  if verbosity >= 2
    println(" ====== Annual average: generalised extreme values.")
    println(" ======    Parameters: μ = ", d_year.μ, "; σ = ", d_year.σ, "; ξ = ", d_year.ξ, ".")
  end

  ## Second with the seasonal variations with respect to the annual averages.
  if verbosity >= 1
    println(" ==== Fitting generalised extreme values models on the seasonal averages.")
  end

  wet_deltas = (wet_avg .- ann_avg) ./ ann_avg
  d_season_wet_params = fit_mle_optim(ExtremeValueDistributions.GeneralizedExtremeValue, wet_deltas, [0.5, 0.5, 0.5], verbose=(verbosity >= 3))
  d_season_wet = Distributions.GeneralizedExtremeValue(d_season_wet_params[1], exp(d_season_wet_params[2]), d_season_wet_params[3]) # ~ \Delta Q_S^W

  if verbosity >= 2
    println(" ====== Wet season average: generalised extreme values.")
    println(" ======    Parameters: μ = ", d_season_wet.μ, "; σ = ", d_season_wet.σ, "; ξ = ", d_season_wet.ξ, ".")
  end

  dry_deltas = (dry_avg .- ann_avg) ./ ann_avg
  d_season_dry_params = fit_mle_optim(ExtremeValueDistributions.GeneralizedExtremeValue, dry_deltas, [0.5, 0.5, 0.5], verbose=(verbosity >= 3))
  d_season_dry = Distributions.GeneralizedExtremeValue(d_season_dry_params[1], exp(d_season_dry_params[2]), d_season_dry_params[3]) # ~ \Delta Q_S^D

  if verbosity >= 2
    println(" ====== Dry season average: generalised extreme values.")
    println(" ======    Parameters: μ = ", d_season_dry.μ, "; σ = ", d_season_dry.σ, "; ξ = ", d_season_dry.ξ, ".")
  end

  ## Finally with the daily variations within each season.
  # Here, working with matrices of values, as each day is considered (and not a single point per year).
  if verbosity >= 1
    println(" ==== Fitting generalised extreme values models on the daily inflow.")
  end

  wet_daily_deltas = vcat([(values(wet_TS[i]) - wet_avg[i]) / wet_avg[i] for i in 1:length(wet_TS)]...)
  d_daily_wet_params = fit_mle_optim(ExtremeValueDistributions.GeneralizedExtremeValue, wet_daily_deltas, [0.5, 0.5, 0.5], verbose=(verbosity >= 3))
  d_daily_wet = Distributions.GeneralizedExtremeValue(d_daily_wet_params[1], exp(d_daily_wet_params[2]), d_daily_wet_params[3]) # ~ \Delta Q_d^W

  if verbosity >= 2
    println(" ====== Daily variations during the wet season: generalised extreme values.")
    println(" ======    Parameters: μ = ", d_daily_wet.μ, "; σ = ", d_daily_wet.σ, "; ξ = ", d_daily_wet.ξ, ".")
  end

  dry_daily_deltas = vcat([(values(dry_TS[i]) - dry_avg[i]) / dry_avg[i] for i in 1:length(dry_TS)]...)
  d_daily_dry_params = fit_mle_optim(ExtremeValueDistributions.GeneralizedExtremeValue, dry_daily_deltas, [0.5, 0.5, 0.5], verbose=(verbosity >= 3))
  d_daily_dry = Distributions.GeneralizedExtremeValue(d_daily_dry_params[1], exp(d_daily_dry_params[2]), d_daily_dry_params[3]) # ~ \Delta Q_d^D

  if verbosity >= 2
    println(" ====== Daily variations during the wet season: generalised extreme values.")
    println(" ======    Parameters: μ = ", d_daily_dry.μ, "; σ = ", d_daily_dry.σ, "; ξ = ", d_daily_dry.ξ, ".")
  end

  ## Done. Prepare the output data structure.
  if verbosity >= 1
    println(" ==== Model fitted. ")
    println(" == ")
  end

  return Adam2015RiverStatisticalModel(d_year, d_season_wet, d_season_dry, d_daily_wet, d_daily_dry, seasonsSeparator)
end

function fit_Adam2015(r::River; kwargs...)
  if haskey(r.ext, "adam2015")
    return r.ext["adam2015"]
  else
    m = fit_Adam2015(r.scenarios; kwargs...)
    r.ext["adam2015"] = m
    return m
  end
end



"""
Prints a report analysing the statistical model with observations made on the raw data.
"""
function printReport(r::River, stat::Adam2015RiverStatisticalModel)
  sc = r.scenarios
  wet = separate(sc, stat.seasons_separator, :Wet)
  dry = separate(sc, stat.seasons_separator, :Dry)
  wet_avg = Float64[mean(values(ts)) for ts in wet]
  dry_avg = Float64[mean(values(ts)) for ts in dry]

  ymean_sc = Float64[mean(values(ts)) for ts in sc]
  wet_smean_sc = (wet_avg .- ymean_sc) ./ ymean_sc
  dry_smean_sc = (dry_avg .- ymean_sc) ./ ymean_sc
  wet_dmean_sc = vcat([(values(wet[i]) .- wet_avg[i]) ./ wet_avg[i] for i in 1:length(wet)]...)
  dry_dmean_sc = vcat([(values(dry[i]) .- dry_avg[i]) ./ dry_avg[i] for i in 1:length(dry)]...)

  # Printing helpers.
  n2s = n -> @sprintf("%7.3f", n)
  function printReport_distr(title::AbstractString, distr::ContinuousUnivariateDistribution, data::Array{Float64, 1})
    # Expected output:
	# 	 ++++ Yearly averages: mean 0.511 vs 0.509
	# 	 ++++              variance 0.014 vs 0.014
	# 	 ++++              kurtosis 0.192 vs 1.689
	#                      ^^^^ (4 characters)
	#         ^^^^^^^^^^^^^^^   length(title), followed by " :"

    title_blank = " " ^ (length(title) - 4)

    println(" ++++ " * title       *     ": mean " * n2s(mean(distr))      * " vs " * n2s(mean(data)))
    println(" ++++ " * title_blank * "  variance " * n2s(var(distr))       * " vs " * n2s(var(data)))
    println(" ++++ " * title_blank * "  skewness " * n2s(skewness(distr))  * " vs " * n2s(skewness(data)))
    println(" ++++ " * title_blank * "  kurtosis " * n2s(kurtosis(distr))  * " vs " * n2s(kurtosis(data)))
  end

  # Validity report.
  println(" ++ Validity report for the statistical model (model vs reality): ")
  printReport_distr("Yearly averages",   stat.year,       ymean_sc)
  printReport_distr("Wet season deltas", stat.season_wet, wet_smean_sc)
  printReport_distr("Dry season deltas", stat.season_dry, dry_smean_sc)
  printReport_distr("Wet daily deltas",  stat.daily_wet,  wet_dmean_sc)
  printReport_distr("Dry daily deltas",  stat.daily_dry,  dry_dmean_sc)

  # Covariance report.
  corr_adj  = c -> (abs(c) < 0.3) ? "weak" : (abs(c) < 0.8) ? "moderate" : "strong"
  corr_sign = c -> (c <= 0.0) ? "negative" : "positive"

  pearson_ymean_wetsmean = cor(ymean_sc, wet_smean_sc)
  pearson_ymean_drysmean = cor(ymean_sc, dry_smean_sc)

  println(" ++ ")
  println(" ++ Covariance report for the statistical model: ")
  println(" ++++ Pearson's coefficient between yearly average and wet season average: " * n2s(pearson_ymean_wetsmean))
  println(" ++++      " * corr_adj(pearson_ymean_wetsmean) * " " * corr_sign(pearson_ymean_wetsmean) * " correlation")
  println(" ++++ Pearson's coefficient between yearly average and dry season average: " * n2s(pearson_ymean_drysmean))
  println(" ++++      " * corr_adj(pearson_ymean_drysmean) * " " * corr_sign(pearson_ymean_drysmean) * " correlation")
end



# Base functions: from the statistical model to actual scenarios.
# Those functions allow scalar *and* vector arguments! (Only useful for Qd, however.)
# 213 days in wet season:     Date(2016, April, 30) - Date(2015, October, 1) + Day(1)
# 152 left in dry season.
adam2015_wet_season_value(QY, QS) = QY * (1 + QS)
adam2015_dry_season_value(QY, QS) = QY * (1 - QS * 213 / 152)
adam2015_wet_day_value(QY, QS, Qd) = QY * (1 + QS) * (1 + Qd)
adam2015_dry_day_value(QY, QS, Qd) = QY * (1 - QS * 213 / 152) * (1 + Qd)



"""
Based on the given statistical `model`, generate one yearly scenarios, with a value each `period`.
This scenario starts at the beginning of the wet season (as per the seasons in `model`), and occur in the year
`currentYear`.
"""
function generateScenarios_Adam2015(model::Adam2015RiverStatisticalModel, period::Period, currentYear::Int; kwargs...)
  # Pick a year average and a season average at random. The wet average has a high likelihood of being negative,
  # so generate until something positive is generated.
  year = rand(getYearlyAverage(model))

  # dry_season = rand(getSeasonAverageDelta(model, :Dry)) # Not fitted, using a fraction of the other season.
  wet_season = rand(getSeasonAverageDelta(model, :Wet))
  while adam2015_wet_season_value(year, wet_season) <= 0.0
    wet_season = rand(getSeasonAverageDelta(model, :Wet))
  end

  # Then work on days.
  bs, es = seasonBoundariesGuarantee(getSeasonSeparator(model), period, currentYear)
  n_wet = length(bs[:Wet] : period : es[:Wet])
  n_dry = length(bs[:Dry] : period : es[:Dry])

  gen_wet = adam2015_wet_day_value(year, wet_season, rand(getSeasonDailyDelta(model, :Wet), n_wet))
  gen_dry = adam2015_dry_day_value(year, wet_season, rand(getSeasonDailyDelta(model, :Dry), n_dry))

  gen = vcat(gen_wet, gen_dry)
  gen[gen .<= 0] = 0.
  return (collect(bs[:Wet] : period : es[:Dry]), gen)
end

generateScenarios_Adam2015(r::River, period::Period, currentYear::Int; kwargs...) =
  generateScenarios_Adam2015(fit_Adam2015(r; kwargs...), period, currentYear)



"""
Computes confidence intervals based on the statistical model.

Assumption: all three stages of probability distributions are not correlated.
This assumption seems mostly valid, according to a few correlation tests.
"""
function computeCI_Adam2015(r::River, model::Adam2015RiverStatisticalModel, level::Float64=.95)
  TS = r.scenarios

  ## Start by computing the actual values to this confidence level using quantiles.
  q_low = (1 - level) / 2
  q_high = 1 - q_low

  # First the annual confidence interval.
  year_low = quantile(model.year, q_low)
  year_mean = mean(model.year)
  year_high = quantile(model.year, q_high)

  # Then the same for the seasonal average normalised deltas.
  wet_season_low = quantile(model.season_wet, q_low)
  wet_season_mean = mean(model.season_wet)
  wet_season_high = quantile(model.season_wet, q_high)

  dry_season_low = quantile(model.season_dry, q_low)
  dry_season_mean = mean(model.season_dry)
  dry_season_high = quantile(model.season_dry, q_high)

  # And finally for the seasonal daily normalised deltas.
  wet_daily_low = quantile(model.daily_wet, q_low)
  wet_daily_mean = mean(model.daily_wet)
  wet_daily_high = quantile(model.daily_wet, q_high)

  dry_daily_low = quantile(model.daily_dry, q_low)
  dry_daily_mean = mean(model.daily_dry)
  dry_daily_high = quantile(model.daily_dry, q_high)

  ## Reconstruct the complete confidence interval as the sum of the previous random variables.
  # Here is used the hypothesis of absence of correlation.
  # Technique:
  #   - start from the annual average bound
  #   - transform one season's distribution based on this to get the seasonal average's bound
  #   - sum those bounds to get the total bound for seasonal variations
  # On average, the daily variations are zero: they are not taken into account when computing the mean.
  # However, these have a larger impact on the bounds of the confidence interval.
  dwet_low  = adam2015_wet_day_value(year_low,  wet_season_low,  wet_daily_low)
  dwet_high = adam2015_wet_day_value(year_high, wet_season_high, wet_daily_high)
  ddry_low  = adam2015_dry_day_value(year_low,  dry_season_low,  dry_daily_low)
  ddry_high = adam2015_dry_day_value(year_high, dry_season_high, dry_daily_high)

  # NOTE: dry_season seems to be replaced by wet_season in Adam's code!?
  dwet_mean = adam2015_wet_season_value(year_mean, wet_season_mean)
  ddry_mean = adam2015_dry_season_value(year_mean, wet_season_mean)

  ## Prepare the outputs, respecting the two seasons.
  wet_begin = begins(model.seasons_separator, :Wet, 1)
  wet_end = ends(model.seasons_separator, :Dry, 1)
  dry_begin = begins(model.seasons_separator, :Dry, 1)
  dry_end = ends(model.seasons_separator, :Dry, 1)
  period = getPeriod(r)

  # The length of the reconstructed series must match that of the time series!
  while length(vcat(wet_begin:period:wet_end, dry_begin:period:dry_end)) > length(timestamp(TS[1]))
    wet_begin += period
  end
  if length(vcat(wet_begin:period:wet_end, dry_begin:period:dry_end)) < length(timestamp(TS[1]))
    error("Reconstructed period too small (" * length(vcat(wet_begin:period:wet_end, dry_begin:period:dry_end)) * ") " *
          "with respect to existing time series (" * length(timestamp(TS[1])) * ").")
  end

  means   = vcat([dwet_mean for d in wet_begin:period:wet_end], [ddry_mean for d in dry_begin:period:dry_end])
  ci_low  = vcat([dwet_low  for d in wet_begin:period:wet_end], [ddry_low  for d in dry_begin:period:dry_end])
  ci_high = vcat([dwet_high for d in wet_begin:period:wet_end], [ddry_high for d in dry_begin:period:dry_end])

  means[means .<= 0.0] = 0.0
  ci_low[ci_low .<= 0.0] = 0.0
  ci_high[ci_high .<= 0.0] = 0.0

  ts_means   = TimeArray(timestamp(TS[1]), means,   colnames(TS[1]), meta(TS[1]))
  ts_ci_low  = TimeArray(timestamp(TS[1]), ci_low,  colnames(TS[1]), meta(TS[1]))
  ts_ci_high = TimeArray(timestamp(TS[1]), ci_high, colnames(TS[1]), meta(TS[1]))

  return ts_means, ts_ci_low, ts_ci_high
end

computeCI_Adam2015(r::River, level::Float64=.95; kwargs...) =
  computeCI_Adam2015(r, fit_Adam2015(r; kwargs...), level)
