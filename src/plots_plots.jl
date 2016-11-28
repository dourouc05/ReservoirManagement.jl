["""
Produce plots about the reservoirs and the solutions.

Available plots:

  * `plotRiverAlongReservoirs`: compares the scenarios for one river for multiple reservoirs.
  * `plotSolutionVsScenarios`: compares a solution to the various scenarios that were used.
  * `plotSolutions`: compares multiple solutions with a reference rule curve.
  * `plotRule`: shows one rule curve (i.e. an array).

Common optional arguments:

  * To control the labels:

    * `title`: title of the plot
    * `lbl_x`: label of the X axis
    * `lbl_y`: label of the Y axis

  * Others:

    * `exportPlotData`: where the function should return only the plot (`false`, default) or the plot and a matrix
      containing all the data used for the plot -- for a CSV export, e.g. -- (`true`).
      The matrix is organised as follows:

      * for the first dimension: first the time-dependent data, then the metadata,
      * for the second dimension: first the solutions, then the scenarios.
"""]



function plotRiverAlongReservoirs(df::DataFrames.DataFrame, period::Period, river::AbstractString;
                                  title::AbstractString="Flow in river " * river,
                                  lbl_x::AbstractString=unitAsString(period), lbl_y::AbstractString="Flow [10^6 m^3/" * unitAsLetter(period) * "]",
                                  exportPlotData::Bool=false)
  plt = plot(df, :x, :y, group=:label, xlabel=lbl_x, ylabel=lbl_y, title=title)
  if ! exportPlotData
    return plt
  else
    ldfg = groupby(df, :label)
    n_groups = length(ldfg)
    n_time_steps = length(ldfg[1][:y])
    dat = Matrix(n_time_steps + 1, n_groups)
    for i in 1:n_groups
      dat[1:n_time_steps, i] = ldfg[i][:y]
      dat[n_time_steps + 1, i] = ldfg[i][:label][1] # All labels are equal within a group (as per groupby).
    end
    return plt, dat
  end
end

function plotRiverAlongReservoirs(df::DataFrames.DataFrame, supp::Dict{AbstractString, Tuple{UnitRange{Int64}, Array{Array{Float64, 1}, 1}}}, period::Period, river::AbstractString;
                                  title::AbstractString="Flow in river " * river,
                                  lbl_x::AbstractString=unitAsString(period), lbl_y::AbstractString="Flow [10^6 m^3/" * unitAsLetter(period) * "]",
                                  exportPlotData::Bool=false)
  labels = sort(vcat(unique(df[:label]), collect(keys(supp))))
  n_lines = length(labels)
  background = RGB{U8}(1.0,1.0,1.0)
  colours = Plots.get_color_palette(:auto, background, n_lines)
  colours_labelled = [labels[i] => colours[i] for i in 1:n_lines]

  # First plot the multiple scenarios, then the other lines (so they are ensured to be on top).
  plt = plot()
  for (label, data) in supp
    # Plot the lines with some transparency.
    raw_col = colours_labelled[label]
    col = RGBA(raw_col.r, raw_col.g, raw_col.b, .15)

    firstLabel = true
    for ys in data[2] # Vector of vectors!
      if firstLabel
        plot!(plt, data[1], ys, line=(col), label=label)
      else
        plot!(plt, data[1], ys, line=(col), label="")
      end
      firstLabel = false
    end
  end
  plot!(plt, df, :x, :y, group=:label, xlabel=lbl_x, ylabel=lbl_y, title=title)

  if ! exportPlotData
    return plt
  else
    ldfg = groupby(vcat(df, supp), :label)
    n_groups = length(ldfg)
    n_time_steps = length(ldfg[1])
    dat = Matrix(n_time_steps + 1, n_groups)
    for i in 1:n_groups
      dat[1:n_time_steps, i] = ldfg[i][:y]
      dat[n_time_steps + 1, i] = ldfg[i][:label][1] # All labels are equal within a group (as per groupby).
    end
    return plt, dat
  end
end

"""
Output a plot for the given reservoirs when studying one particular river.
"""
function plotRiverAlongReservoirs(reservoirs::Array{Reservoir, 1}, period::Period, timePeriod::UnitRange{Int64}, river::AbstractString; kwargs...)
  df, supp = getDataFrameForRiver(reservoirs, timePeriod, river, multipleScenarios=:FullSplit)
  # plotRiverAlongReservoirs(df, supp, period, river; kwargs...)
  plotRiverAlongReservoirs(df, period, river; kwargs...)
end

"""
Plots a solution with the solution to each of the scenarios that was used for the optimisation (this only makes sense
for `:BasicRuleCurveSolver`), with a comparison to the current rule curve `currentTarget` (whose length is used to
select the relevant parts of the individual scenarios' solutions).
These scenarios may be coloured according to some characteristics with the `colouring` parameter.

Available `colouring` techniques:

  * `:YearlyAverage`: average over the whole input scenarios
  * `:WetSeasonAverage`: average over the wet season of each scenario (seasons are given by `seasonsSeparator`)
  * `:DrySeasonAverage`: average over the dry season of each scenario (seasons are given by `seasonsSeparator`)
  * `:DriestMonth`: average over the driest month within each scenario
  * `:DriestThreeMonths`: average over the three driest months within each scenario
  * `:DriestSixMonths`: average over the six driest months within each scenario

If `exportPlotData` is `true`, then the function returns two things: the plot and a matrix of the data that are plotted.
Each line of the matrix contains the information for one line in the plot (be it a solution or a scenario). The first
`countTimeSteps(sol)` columns contain the value at each time step. The next columns are:

  * for solutions: their name. A zero-length string indicates a scenario.
  * for scenarios: the order for the given `colouring` order. `1` is the lowest, and indices increase. `-1` indicates a
    solution.
  * for scenarios: the average that is used for the `colouring` as a real number. `-1` indicates a solution.
  * for scenarios: `0` for support scenarios, `1` for nonsupport scenarios. `-1` indicates a solution.

"""
function plotSolutionVsScenarios(sol::ReservoirSolution, currentTarget::Array{Float64, 1};
                               colouring::Symbol=:YearlyAverage, seasonsSeparator::SeasonsSeparator=HydrologicalWesternEurope,
                               title::AbstractString="Reservoir storage for solver " * sol.solverDescription * ": \ncomparing scenarios",
                               lbl_x::AbstractString=unitAsString(getPeriod(sol)), lbl_y::AbstractString="Reservoir storage [m^3]",
                               lbl_ct::AbstractString="Current rule curve", lbl_ns::AbstractString="Optimised rule curve",
                               col_ct::Colorant=colorant"purple", col_ns::Colorant=colorant"steel blue",
                               exportPlotData::Bool=false)
  ## Basic definitions.
  # Define the colouring algorithms.
  transform = Dict(
    :YearlyAverage     => x -> Float64[mean(values(x[sc])) for sc in 1:size(x, 1)],
    :WetSeasonAverage  => x -> Float64[mean(values(collect(x[sc], seasonsSeparator, :Wet))) for sc in 1:size(x, 1)],
    :DrySeasonAverage  => x -> Float64[mean(values(collect(x[sc], seasonsSeparator, :Dry))) for sc in 1:size(x, 1)],
    :DriestMonth       => x -> Float64[mean(sort(values(collapse(x[sc], month, first, sum)))[1]) for sc in 1:size(x, 1)],
    :DriestThreeMonths => x -> Float64[mean(sort(values(collapse(x[sc], month, first, sum)))[1:3]) for sc in 1:size(x, 1)],
    :DriestSixMonths   => x -> Float64[mean(sort(values(collapse(x[sc], month, first, sum)))[1:6]) for sc in 1:size(x, 1)],
  )
  colourings = collect(keys(transform))

  # Check the given arguments are within constraints. Only :BasicRuleCurveSolver is supported: other solvers
  # do not use scenarios, at least not in a way that is meaningful for this plot (safeMinimumRuleCurve considers
  # a series of scenarios for each time step in output, that would make little sense here).
  @assert sol.solver == :BasicRuleCurveSolver "Wrong kind of solver! This function only accepts :BasicRuleCurveSolver."
  @assert in(colouring, colourings) "Unrecognised colouring."

  ## Pretreatment before plotting.
  # Prepare the data for plotting.
  xs = 1:length(currentTarget)
  ts = vcat([sol.solution.levels[i, xs] for i in 1:size(sol.solution.levels, 1)]...)
  ts_colours = colormap("RdBu", size(ts, 1)) # From red to blue.

  # According to the colouring criterion, sort the scenarios when giving them a colour.
  ins = totalInflow(sol.reservoir, 1:countScenarios(sol.reservoir), asTimeSeries=true)
  avgs = transform[colouring](ins)
  ts_order = sortperm(avgs, rev=true) # Lowest values at the end!
  choose_colour = (index::Int) -> ts_colours[ts_order[index]]

  # Determine which scenarios are support vectors, i.e. the total solution touches the solution to this scenario.
  support = falses(countScenarios(sol.reservoir))
  support_count_threshold = 3 # Number of points that must touch the solution to be considered as support scenario. TODO: parameter to the function. Default to 1?
  for s in 1:countScenarios(sol.reservoir)
    if sum(sol.solution.levels[s, xs] .== sol.solution.solution[xs]) >= support_count_threshold
      #                               ^^^
      # Can do a crude equality test: the values in the solution are copied from those of scenarios (no extra computation).
      # No need for a difference and a threshold here.
      support[s] = true
    end
  end

  ## Actual plotting.
  # Plot the different solutions.
  plt = plot(xs, currentTarget[xs],         label=lbl_ct, linewidth=3, linecolor=col_ct, xlabel=lbl_x, ylabel=lbl_y, title=title)
  plot!(plt, xs, sol.solution.solution[xs], label=lbl_ns, linewidth=5, linecolor=col_ns)

  # Plot the current target and the optimised target.
  for s in 1:countScenarios(sol.reservoir)
    col = choose_colour(s)
    lbl = (col == ts_colours[1]) ? "driest" : (col == ts_colours[end]) ? "wettest" : ""
    plot!(plt, xs, vec(sol.solution.levels[s, xs]), label=lbl, line=(choose_colour(s), support[s] ? :solid : :dash))
  end

  # Return what you should: plot plus sometimes data.
  if ! exportPlotData
    return plt
  else
    dat = Matrix(length(xs) + 4, 2 + size(sol.solution.levels, 1))
    c_name = length(xs) + 1
    c_order = length(xs) + 2
    c_avg = length(xs) + 3
    c_support = length(xs) + 4

    # First deal with solutions.
    dat[xs, 1] = currentTarget[xs]
    dat[c_name, 1] = lbl_ct
    dat[xs, 2] = sol.solution.solution[xs]
    dat[c_name, 2] = lbl_ns
    dat[c_order:c_support, 1:2] = -1

    # Second deal with scenarios.
    dat[xs, 3:2 + size(sol.solution.levels, 1)] = sol.solution.levels[:, xs]'
    dat[c_name, 3:2 + size(sol.solution.levels, 1)] = -1
    dat[c_order, 3:2 + size(sol.solution.levels, 1)] = ts_order
    dat[c_avg, 3:2 + size(sol.solution.levels, 1)] = avgs
    dat[c_support, 3:2 + size(sol.solution.levels, 1)] = support

    return plt, dat
  end
end

function plotSolutions(df::DataFrames.DataFrame, period::Period;
              title::AbstractString="Reservoir storage depending on optimisation method",
              lbl_x::AbstractString=unitAsString(period), lbl_y::AbstractString="Reservoir storage [m^3]",
              exportPlotData::Bool=false)
  plt = plot(df, :x, :y, group=:label, xlabel=lbl_x, ylabel=lbl_y, title=title)

  if ! exportPlotData
    return plt
  else
    ldfg = groupby(df, :label)
    n_groups = length(ldfg)
    n_time_steps = size(ldfg[1], 1)
    dat = Matrix(n_time_steps + 1, n_groups)
    for i in 1:n_groups
      dat[1:n_time_steps, i] = ldfg[i][:y]
      dat[n_time_steps + 1, i] = ldfg[i][:label][1] # All labels are equal within a group (as per groupby).
    end
    return plt, dat
  end
end

function plotSolutions(solutions::Array{ReservoirSolution, 1}, currentTarget::Array{Float64, 1}; period::Union{Void, Period}=nothing, kwargs...)
  @assert length(solutions) != 0 "No solutions to plot. Please provide at least one solution."
  if period == nothing
    period = getPeriod(solutions[1])
  end

  xs = 1:length(currentTarget)
  df = vcat(getDataFrameForSolutions(solutions, xs), DataFrame(x=xs, y=currentTarget, label="Current rule curve"))
  return plotSolutions(df, period; kwargs...)
end

plotRule(currentTarget::Array{Float64, 1}, period::Period; kwargs...) =
  plotSolutions(DataFrame(x=1:length(currentTarget), y=currentTarget, label="Current rule curve"), period; kwargs...)
