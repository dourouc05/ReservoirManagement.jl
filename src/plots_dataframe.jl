["""
Defines generic plotting facilities.
"""]



"""
Produces a `DataFrame` for the given reservoirs when studying one particular river. Multiple modes are possible to handle the various scenarios (parameter `multipleScenarios`):
  * `:None` does not allow multiple scenarios for any given river, and will return an error if this hypothesis is proven wrong.
  * `:Ignore` only considers the first scenario for any given river.
  * `:Merge` includes all scenarios inside the resulting data frame, with the index of the scenario appended in the column name.
  * `:Split` separates the first scenarios and then the other ones for each reservoir (e.g., to plot the first output with the legend, and add the others).
  * `:FullSplit` separates the reservoirs with one scenario from those having more than one.
    The second output links the time steps (only integer values) and a vector of vector of values for each time step; these tuples are indexed by .
    The return values has then the type: `Tuple{DataFrames.DataFrame, Dict{AbstractString, Tuple{UnitRange{Int64}, Array{Array{Float64, 1}, 1}}}}`.
"""
function getDataFrameForRiver(reservoirs::Array{Reservoir, 1}, timePeriod::UnitRange{Int64}, river::AbstractString; multipleScenarios::Symbol=:None)
  for r in reservoirs
    if ! hasRiver(r, river)
      error("Unknow river " * river * " for reservoir " * string(r) * ".")
    end
  end

  if multipleScenarios == :None
    for i in 1:length(reservoirs)
      if countScenarios(get(getRiver(reservoirs[i], river))) > 1
        error("At least one reservoir (id: " * string(i) * ") has more than one scenario for the river " * river)
      end
    end
    return getDataFrameForRiver(reservoirs, timePeriod, river, multipleScenarios=:Ignore)
  elseif multipleScenarios == :Ignore
    return vcat([DataFrame(x=timePeriod, y=values(get(getRiver(reservoirs[i], river)).scenarios[1][timePeriod]), label=reservoirs[i].variant) for i in 1:length(reservoirs)]...)
  elseif multipleScenarios == :Merge
    label = (i, j) -> reservoirs[i].variant * " (" * string(j) * ")"
    return vcat([
        let r = get(getRiver(reservoirs[i], river))
          vcat([
            DataFrame(x=timePeriod, y=values(r.scenarios[j][timePeriod]), label=label(i, j))
            for j in 1:countScenarios(r)
          ]...)
        end
        for i in 1:length(reservoirs)
      ]...)
  elseif multipleScenarios == :Split
    main = getDataFrameForRiver(reservoirs, timePeriod, river, multipleScenarios=:Ignore)
    copies = Dict{AbstractString, Tuple{UnitRange{Int64}, Array{Array{Float64, 1}, 1}}}()
    for i in 1:length(reservoirs)
      r = get(getRiver(reservoirs[i], river, nullable=false))
      if countScenarios(r) > 1
        copies[reservoirs[i].variant] = (timePeriod, Array{Float64, 1}[values(r.scenarios[j][timePeriod]) for j in 2:countScenarios(r)])
      end
    end
    return main, copies
  elseif multipleScenarios == :FullSplit
    main = DataFrame()
    copies = Dict{AbstractString, Tuple{UnitRange{Int64}, Array{Array{Float64, 1}, 1}}}()
    for i in 1:length(reservoirs)
      r = get(getRiver(reservoirs[i], river))
      if countScenarios(r) == 1
        main = vcat(main, DataFrame(x=timePeriod, y=values(r.scenarios[1][timePeriod]), label=reservoirs[i].variant))
      else
        copies[reservoirs[i].variant] = (timePeriod, Array{Float64, 1}[values(r.scenarios[j][timePeriod]) for j in 1:countScenarios(r)])
      end
    end
    return main, copies
  else
    error("Unknown handling requirement.")
  end
end

"""
Produces a `DataFrame` for the given solutions, encoding the solution for each algorithm.
"""
getDataFrameForSolutions(solutions::Array{ReservoirSolution, 1}, timePeriod::UnitRange{Int64}) =
  vcat(DataFrame[DataFrame(x=timePeriod, y=solutions[i].solution.solution[timePeriod], label=getNameOrElse(solutions[i])) for i in 1:length(solutions)]...)



unitAsLetter(period::Period) =
  if period == Day(1)
    "d"
  elseif period == Week(1)
    "w"
  elseif period == Month(1)
    "m"
  elseif period == Year(1)
    "y"
  else
    error("Unrecognised period: ", period)
  end

unitAsString(period::Period) =
  if period == Day(1)
    "day"
  elseif period == Week(1)
    "week"
  elseif period == Month(1)
    "month"
  elseif period == Year(1)
    "year"
  else
    error("Unrecognised period: ", period)
  end
