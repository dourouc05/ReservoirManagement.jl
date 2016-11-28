["""
Solves the inverse problem to find the minimum inflow so the reservoir can meet its functions.
"""]



immutable FunctionsInverseOptions <: InverseOptions
  objective::Symbol
end

immutable FunctionsInverseSolution <: InverseSolution
  solution::Array{Float64, 1}
  total::Float64
end

immutable FunctionsInverseSolverStatistics <: InverseSolverStatistics
  total::Float64
  modelling::Float64
  writingModelFile::Float64
  solving::Float64
  postprocessing::Float64
end



"""
Compute the worst case for a given solution: if the input is lower than the output of this function, then
the input solution cannot be feasible.

Inputs:

  * `solution`: a solution to evaluate for feasibility (implicitly, a reservoir)
  * `objective`: a choice of objective function.
    * `:Total`: minimise the total precipitations
  * `verbosity`: level of output.
    * `-1`: no output
    * `0`: no output, except for solver's warnings
    * `1`: summary output
    * `2`: detailed output of this algorithm's details
    * `3`: detailed output of this algorithm's and its sub-algorithms' details
    * `4`: including solver output
"""
function worstCaseFunctions(solution::ReservoirSolution;
                            objective::Symbol=:Total, useOutputs::Array{Symbol, 1}=Symbol[],
                            name="", output="", verbosity::Int=1)
  reservoir = solution.reservoir

  if objective != :Total
    error("Unrecognised objective.")
  end

  if verbosity >= 1
    println(" == ")
    println(" ==== Optimisation model for inflow: ", getVariant(reservoir))
  end

  time_start = float(time_ns())

  # Prepare the data.
  time_horizon = length(solution.solution.solution)

  if verbosity >= 2
    println(" ====== Time horizon: ", time_horizon)
  end

  # Create the optimisation model.
  m = Model(solver=GurobiSolver(OutputFlag=verbosity >= 4))
  @variable(m, getMinimumCapacity(reservoir) <= storage[0:time_horizon] <= getMaximumCapacity(reservoir)) # 10^6 m^3
  @variable(m, 0 <= inflow[1:time_horizon] <= getMaximumCapacity(reservoir)) # 10^6 m^3
  @variable(m, release[1:time_horizon] >= 0) # 10^6 m^3/time step

  @objective(m, Min, sum(inflow))

  for t = 1:time_horizon
    @constraint(m, storage[t] == storage[t - 1] - hlp_JuMP_purposes_discharge(m, reservoir, getPeriod(reservoir)) + inflow[t] - release[t])
    @constraint(m, storage[t] == solution.solution.solution[t])
    hlp_JuMP_outputs_constraint(m, release[t], storage[t], reservoir, useOutputs)
  end

  time_modelled = float(time_ns())

  # Services around the model.
  if output != ""
    writeLP(m, output)
  end
  time_output = float(time_ns())

  status = solve(m); time_solved = float(time_ns())
  if status != :Infeasible
    inflow = getvalue(inflow)
    totalInflow = getobjectivevalue(m)
    feasible = true

    if verbosity >= 1
      println(" ==== Total inflow: ", totalInflow * 1000000, " m^3.")
    end
  else
    inflow = zeros(time_horizon) / 0. # NaN as 0/0.
    totalInflow = inflow[1]
    feasible = false

    if verbosity >= 1
      println(" ==== Infeasible!") # What the heck? Where did you get that solution from?
    end
  end

  time_done = float(time_ns())

  if verbosity >= 1
    println(" ==== Total time: ", (time_solved - time_start) / 1000000, "ms")
    if verbosity >= 2
      println(" ====== Time to create the JuMP model: ", (time_modelled - time_start) / 1000000, "ms")
      println(" ====== Time to write the LP model if required: ", (time_output - time_modelled) / 1000000, "ms")
      println(" ====== Time to actually solve the model: ", (time_solved - time_output) / 1000000, "ms")
      println(" ====== Time to prepare the return values: ", (time_done - time_solved) / 1000000, "ms")
    end
    println(" == ")
  end

  options = FunctionsInverseOptions(objective)
  solution = FunctionsInverseSolution(inflow, totalInflow)
  time = FunctionsInverseSolverStatistics((time_solved - time_start) / 1000000, (time_modelled - time_start) / 1000000, (time_output - time_modelled) / 1000000, (time_solved - time_output) / 1000000, (time_done - time_solved) / 1000000)
  solName = (name == "") ? Nullable{AbstractString}() : Nullable(name)
  return ReservoirSolution(reservoir, feasible, solution, solName, :FunctionsInverseSolver, "inverse inflow from rule curve", options, time)
end
