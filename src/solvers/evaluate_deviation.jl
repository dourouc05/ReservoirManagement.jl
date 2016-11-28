["""
Computes the maximum deviation from some reference scenario a solution can withstand.
"""]



immutable DeviationOptions <: EvaluationOptions
  useOutputs::Array{Symbol, 1}
end

immutable DeviationSolution <: EvaluationSolution
  deviation::Float64
end

immutable DeviationSolverStatistics <: EvaluationSolverStatistics
  total::Float64
  modelling::Float64
  solving::Float64
  writingModelFile::Float64
end



"""
Computes the maximum deviation from some reference inflow scenario a solution can withstand.
For example, the input scenario may be the result of the inverse model on some preexisting solution.

  * `solution`: a solution to evaluate for feasibility (implicitly, a reservoir)
  * `inflow`: the reference inflow

The available ways of releasing water are configurable. If all of them are disabled (by default), there is no bound
on the output. Otherwise, each dam output present in the vector `useOutputs` is enabled. The outputs are identified
with their id.

Finally, `verbosity` controls the level of output from this solver:

  * `-1`: no output
  * `0`: no output, except for solver's warnings
  * `1`: summary output
  * `2`: detailed output of this algorithm's details
  * `3`: detailed output of this algorithm's and its sub-algorithms' details
  * `4`: including solver output
"""
function maximumDeviation(solution::ReservoirSolution, inflow::Array{Float64, 1}, name::AbstractString="";
                          useOutputs::Array{Symbol, 1}=Symbol[],
                          output="", verbosity::Int=1)
  reservoir = getReservoir(solution)
  if verbosity >= 1
    println(" == ")
    println(" ==== Maximum deviation analysis: ", getVariant(reservoir))
  end

  time_start = float(time_ns())

  ## Deal with input data.
  time_horizon = length(inflow)
  minLevel = getSolution(solution)
  period = getPeriod(reservoir)

  if verbosity >= 2
    println(" ====== Input time horizon: ", time_horizon)
  end

  ## Create the optimisation model.
  m = Model(solver=GurobiSolver(OutputFlag=(verbosity >= 4)))
  @variable(m, delta >= 0) # 10^6 m^3
  @variable(m, getMinimumCapacity(reservoir) <= storage[0:time_horizon] <= getMaximumCapacity(reservoir)) # 10^6 m^3
  @variable(m, release[1:time_horizon] >= 0) # 10^6 m^3/time step

  # Objective function.
  @objective(m, Max, delta)

  # Main set of constraints.
  for t = 1:time_horizon
    # Implement the state equation.
    @constraint(m, storage[t] == storage[t - 1] + inflow[t]
                                 - hlp_JuMP_purposes_discharge(m, reservoir, period) - release[t] - delta)
    hlp_JuMP_outputs_constraint(m, release[t], storage[t], reservoir, useOutputs)

    # Implement the rule curve.
    @constraint(m, storage[t] >= minLevel[t])
  end

  time_modelled = float(time_ns())

  ## Services around the model.
  if output != ""
    writeLP(m, output)
  end
  time_output = float(time_ns())

  status = solve(m, suppress_warnings=verbosity < 0)
  if status != :Infeasible
    max_deviation = getvalue(delta)
    feasible = true

    if verbosity >= 1
      println(" ==== Maximum deviation: ", getobjectivevalue(m) * 1000000, " m^3 over ", time_horizon * getPeriod(reservoir), ".")
    end
  else
    if verbosity >= 1
      println(" ==== Infeasible!")
    end
    feasible = false
  end

  time_done = float(time_ns())

  # Compute the timings (in ms).
  total_time = (time_done - time_start) / 1000000
  modelling_time = (time_modelled - time_start) / 1000000
  writing_time = (time_output - time_modelled) / 1000000
  solving_time = (time_done - time_output) / 1000000

  # Final output.
  if verbosity >= 1
    println(" ==== Total time: ", total_time, "ms")
    if verbosity >= 2
      println(" ====== Time to create the JuMP model: ", modelling_time, "ms")
      println(" ====== Time to write the LP model if required: ", writing_time, "ms")
      println(" ====== Time to actually solve the model: ", solving_time, "ms")
    end
    println(" == ")
  end

  # Prepare the data structures for the outside world.
  options = DeviationOptions(useOutputs)
  solution = DeviationSolution(max_deviation)
  time = DeviationSolverStatistics(total_time, modelling_time, writing_time, writing_time)
  solName = (name == "") ? Nullable{AbstractString}() : Nullable(name)
  return ReservoirSolution(reservoir, feasible, solution, solName, :DeviationSolver, "maximum deviation from inflow", options, time)
end
