["""
Evaluates the shortage or surplus of water when trying to meet all purposes. To determine the reason a solver fails,
this might have two outputs: either there is not enough water (impossible to meet purposes), or too much (flood control
issues).
"""]



immutable ShortageOptions <: EvaluationOptions
  useOutputs::Array{Symbol, 1}
  forceInitialCondition::Bool
end

immutable ShortageSolution <: EvaluationSolution
  infeasibilityCause::Symbol # Either :Shortage or :Surplus, or :MostlyShortage or :MostlySurplus,
  # or :None if the solver found no infeasibility cause (i.e. the given rule curve is feasible).
  delta::Array{Float64, 1} # Delta over time: shortage if positive, surplus if negative.
end

immutable ShortageSolverStatistics <: EvaluationSolverStatistics
  total::Float64
  modelling::Array{Float64, 1}
  writingModelFile::Array{Float64, 1}
  solving::Array{Float64, 1}
end



"""
Computes the shortage index based on the result of a `purposeShortage` analysis (given in `s`),
for the given `purposes`. This index compares the missing water to meet all purposes with the total need for water,
year by year:

    SI = (1 / n_years) * \sum_{i \in years} (shortage_i / total_needs_i)^2
    shortage_i = demand_i - supplied_i

This index makes sense only for deterministic water withdrawal: otherwise, it is impossible to characterise the total
water that is needed for a given set of purposes. Also, when multiple purposes are considered for the optimisation,
it does not make sense to compute the index with respect to only one, as there is no way to allot the missing water
between purposes. (This would need some prioritisation between the purposes, which is not available.)

Reference: http://www.hec.usace.army.mil/publications/IHDVolumes/IHD-8.pdf, page 66.
"""
function computePurposeShortageIndex(s::ReservoirSolution, purposes::Array{Purpose, 1})
  # Check whether this function is applicable: only makes sense for deterministic water withdrawal purposes.
  for p in purposes
    if ! isWaterWithdrawal(p) || ! isDeterministic(p)
      error("A purpose is not fully determined by a deterministic water withdrawal: " *
            string(getId(p)) * " (type: " * string(typeof(p)) * ")")
    end
  end

  if ! (typeof(getSolutionObject(s)) <: ShortageSolution)
    error("Not the solution of a shortage solver, cannot compute the index. Use purposeShortage().")
  end

  solution = getSolutionObject(s)
  reservoir = getReservoir(s)
  period = getPeriod(reservoir)

  # If there is no shortage, return 0.
  if solution.infeasibilityCause != :Shortage
    return 0.
  end

  # Determine the ranges of values to aggregate so that the index computations are performed year by year.
  time_horizon = length(solution.delta)
  shortage = solution.delta
  shortage[shortage .<= 0.] = 0.
  demands = [sum([getNeed(p, period) for p in getPurposes(reservoir)]) for t in 1:time_horizon]

  dates = getTimePoints(reservoir)
  if length(dates) > time_horizon
    dates = dates[1:time_horizon]
  end
  index_beginning = [1] # The first year begins at index 1; then compute the beginning of each following year.
  for i in eachindex(dates)
    if i == 1
      continue
    end

    if dates[i] >= dates[last(index_beginning)] + Year(1)
      push!(index_beginning, i)
    end
  end
  if length(index_beginning) == 1 # Not a full year: consider all values as a part of this "year".
    push!(index_beginning, length(dates) + 1)
  end

  # Actually compute the index, making one value per year.
  shortages_yearly = Float64[sum([shortage[j] for j in index_beginning[i] : index_beginning[i + 1] - 1]) for i in 1:length(index_beginning)-1]
  needs_yearly = Float64[sum([demands[j] for j in index_beginning[i] : index_beginning[i + 1] - 1]) for i in 1:length(index_beginning)-1]

  return mean([(shortages_yearly[i] / needs_yearly[i])^2 for i in eachindex(shortages_yearly)])
end
computePurposeShortageIndex(s::ReservoirSolution, p::Purpose) = computePurposeShortageIndex(s, [p])

"""
Computes the shortage vulnerability index based on the result of a `purposeShortage` analysis (given in `s`),
for the given `purposes`. This index compares the missing water to meet all purposes with the total need for water
continuously:

    VI = (1 / n_time_steps) * \sum_t (shortage_t / total_needs_t)^2
    shortage_t = demand_t - supplied_t

This index makes sense only for deterministic water withdrawal: otherwise, it is impossible to characterise the total
water that is needed for a given set of purposes. Also, when multiple purposes are considered for the optimisation,
it does not make sense to compute the index with respect to only one, as there is no way to allot the missing water
between purposes. (This would need some prioritisation between the purposes, which is not available.)

The difference with the shortage index (`computePurposeShortageIndex`) is that the quotient is evaluated at each
time step, and not only once per year with the total shortage and needs for this year.
"""
function computePurposeVulnerabilityIndex(s::ReservoirSolution, purposes::Array{Purpose, 1})
  # Check whether this function is applicable: only makes sense for deterministic water withdrawal purposes.
  for p in purposes
    if ! isWaterWithdrawal(p) || ! isDeterministic(p)
      error("A purpose is not fully determined by a deterministic water withdrawal: " *
            string(getId(p)) * " (type: " * string(typeof(p)) * ")")
    end
  end

  if ! (typeof(getSolutionObject(s)) <: ShortageSolution)
    error("Not the solution of a shortage solver, cannot compute the index. Use purposeShortage().")
  end

  solution = getSolutionObject(s)
  reservoir = getReservoir(s)
  period = getPeriod(reservoir)

  # If there is no shortage, return 0.
  if solution.infeasibilityCause != :Shortage
    return 0.
  end

  # Actually compute the index.
  time_horizon = length(solution.delta)
  shortage = solution.delta
  shortage[shortage .<= 0.] = 0.
  demands = [sum([getNeed(p, period) for p in getPurposes(reservoir)]) for t in 1:time_horizon]

  return mean([(shortage[t] / demands[t])^2 for t in 1:time_horizon])
end
computePurposeVulnerabilityIndex(s::ReservoirSolution, p::Purpose) = computePurposeVulnerabilityIndex(s, [p])


## First model: check for shortage of water.
function __modelShortage(reservoir::Reservoir, time_horizon::Int, forceInitialCondition::Bool,
                         inflow::Array{Float64, 1}, period::Period, useOutputs::Array{Symbol, 1},
                         minLevel::Array{Float64, 1},
                         verbosity::Int)
  time_start = float(time_ns())

  m = Model(solver=GurobiSolver(OutputFlag=(verbosity >= 4)))
  @variable(m, delta[1:time_horizon] >= 0) # 10^6 m^3
  @variable(m, getMinimumCapacity(reservoir) <= storage[0:time_horizon] <= getMaximumCapacity(reservoir)) # 10^6 m^3
  @variable(m, release[1:time_horizon] >= 0) # 10^6 m^3/time step

  # Objective function.
  @objective(m, Min, sum(delta)) # delta >= 0

  # Initial condition: minimum reservoir level (shortage).
  if forceInitialCondition
    @constraint(m, storage[0] == getMinimumCapacity(reservoir))
  end

  # Main set of constraints.
  for t = 1:time_horizon
    # Implement the state equation.
    @constraint(m, storage[t] == storage[t - 1] + inflow[t] - hlp_JuMP_purposes_discharge(m, reservoir, period)
                                 - release[t] + delta[t])
    hlp_JuMP_outputs_constraint(m, release[t], storage[t], reservoir, useOutputs)

    # Implement the rule curve.
    @constraint(m, storage[t] >= minLevel[t])

    # At any time step, cannot lack more water than what is needed. (Otherwise, the solver might group all lacks
    # at one time step, and only the objective value then makes sense.) Respecting the rule curve might require
    # some water.
    delta_rule_curve = minLevel[t] - (t == 1 ? getMinimumCapacity(reservoir) : minLevel[t - 1])
    delta_rule_curve_plus = delta_rule_curve >= 0. ? delta_rule_curve : 0.
    purposes_discharge = hlp_JuMP_purposes_discharge(m, reservoir, period)
    @constraint(m, delta[t] <= hlp_JuMP_purposes_discharge(m, reservoir, period) + delta_rule_curve_plus)
  end

  time_modelled = float(time_ns())

  # Output the model.
  # if output != ""
  #   writeLP(m, output)
  # end
  time_output = float(time_ns())

  # Solve the model. It might be infeasible if the problem does not exclusively come from a water shortage
  # (i.e. it is due to either exclusively a surplus, or a combination of both).
  status = solve(m, suppress_warnings=verbosity < 0)
  time_solved = float(time_ns())

  objective = (status == :Optimal) ? getobjectivevalue(m) : 0. / 0.
  deltas = (status == :Optimal) ? getvalue(delta) : zeros(time_horizon) / 0.

  return (status, objective, deltas, (time_modelled - time_start) / 1000000, (time_output - time_modelled) / 1000000, (time_solved - time_output) / 1000000)
end

## Second model: check for surplus of water.
function __modelSurplus(reservoir::Reservoir, time_horizon::Int, forceInitialCondition::Bool,
                        inflow::Array{Float64, 1}, period::Period, useOutputs::Array{Symbol, 1},
                        minLevel::Array{Float64, 1},
                        verbosity::Int)
  time_start = float(time_ns())

  m = Model(solver=GurobiSolver(OutputFlag=(verbosity >= 4)))
  @variable(m, delta[1:time_horizon] <= 0) # 10^6 m^3
  @variable(m, getMinimumCapacity(reservoir) <= storage[0:time_horizon] <= getMaximumCapacity(reservoir)) # 10^6 m^3
  @variable(m, release[1:time_horizon] >= 0) # 10^6 m^3/time step

  # Objective function.
  @objective(m, Max, sum(delta)) # delta <= 0

  # Initial condition: maximum reservoir level for surplus.
  if forceInitialCondition
    @constraint(m, storage[0] == getMaximumCapacity(reservoir))
  end

  # Main set of constraints.
  for t = 1:time_horizon
    # Implement the state equation.
    @constraint(m, c, storage[t] == storage[t - 1] + inflow[t] - hlp_JuMP_purposes_discharge(m, reservoir, period)
                                  - release[t] + delta[t])
    hlp_JuMP_outputs_constraint(m, release[t], storage[t], reservoir, useOutputs)

    # Implement the rule curve.
    @constraint(m, storage[t] >= minLevel[t])

    # No need for a bound the on deltas: they are already negative, and only withdraw water from the reservoir
    # (i.e. the purposes are not endangered).
  end

  time_modelled = float(time_ns())

  # Output the model.
  # if output != ""
  #   writeLP(m, output)
  # end
  time_output = float(time_ns())

  # Solve the model.
  status = solve(m, suppress_warnings=verbosity < 0)
  time_solved = float(time_ns())

  objective = (status == :Optimal) ? getobjectivevalue(m) : 0. / 0.
  deltas = (status == :Optimal) ? getvalue(delta) : zeros(time_horizon) / 0.

  return (status, objective, deltas, (time_modelled - time_start) / 1000000, (time_output - time_modelled) / 1000000, (time_solved - time_output) / 1000000)
end

## Third model: look for more interesting cases, which cannot be explained by the two previous models.
function __modelMixed(reservoir::Reservoir, time_horizon::Int, forceInitialCondition::Bool,
                      inflow::Array{Float64, 1}, period::Period, useOutputs::Array{Symbol, 1},
                      minLevel::Array{Float64, 1},
                      verbosity::Int)
  time_start = float(time_ns())

  m = Model(solver=GurobiSolver(OutputFlag=(verbosity >= 4)))
  @variable(m, delta[1:time_horizon]) # 10^6 m^3
  @variable(m, deltaAbs[1:time_horizon] >= 0) # 10^6 m^3; absolue value of the previous
  @variable(m, getMinimumCapacity(reservoir) <= storage[0:time_horizon] <= getMaximumCapacity(reservoir)) # 10^6 m^3
  @variable(m, release[1:time_horizon] >= 0) # 10^6 m^3/time step

  # Objective function.
  @objective(m, Min, sum(delta)) # delta >= 0

  # Initial condition: average between minimum and maximum reservoir level (shortage and surplus).
  if forceInitialCondition
    @constraint(m, 2 * storage[0] == getMinimumCapacity(reservoir) + getMaximumCapacity(reservoir))
  end

  # Main set of constraints.
  for t = 1:time_horizon
    # Implement the state equation.
    @constraint(m, storage[t] == storage[t - 1] + inflow[t] - hlp_JuMP_purposes_discharge(m, reservoir, period)
                                - release[t] + delta[t])
    hlp_JuMP_outputs_constraint(m, release[t], storage[t], reservoir, useOutputs)

    # Implement the rule curve.
    @constraint(m, storage[t] >= minLevel[t])

    # At any time step, cannot lack more water than what is needed. (Otherwise, the solver might group all lacks
    # at one time step, and only the objective value then makes sense.) Respecting the rule curve might require
    # some water.
    # However, no lower bound, let it at -\infty: for water surplus, might need to remove more water than this.
    delta_rule_curve = minLevel[t] - (t == 1 ? getMinimumCapacity(reservoir) : minLevel[t - 1])
    delta_rule_curve_plus = delta_rule_curve >= 0. ? delta_rule_curve : 0.
    @constraint(m, delta[t] <= hlp_JuMP_purposes_discharge(m, reservoir, period) + delta_rule_curve_plus)

    # Implement the absolute value for the objective function.
    @constraint(m, deltaAbs[t] >= delta[t])
    @constraint(m, deltaAbs[t] >= - delta[t])
  end

  time_modelled = float(time_ns())

  # Output the model.
  # if output != ""
  #   writeLP(m, output)
  # end
  time_output = float(time_ns())

  # Solve the model. It might be infeasible if the problem does not exclusively come from a water shortage
  # (i.e. it is due to either exclusively a surplus, or a combination of both).
  status = solve(m, suppress_warnings=verbosity < 0)
  time_solved = float(time_ns())

  objective = (status == :Optimal) ? getobjectivevalue(m) : 0. / 0.
  deltas = (status == :Optimal) ? getvalue(delta) : zeros(time_horizon) / 0.

  return (status, objective, deltas, (time_modelled - time_start) / 1000000, (time_output - time_modelled) / 1000000, (time_solved - time_output) / 1000000)
end



"""
Evaluates the shortage of water when trying to meet all purposes, i.e. the quantity of supplementary water (with
respect to the inflow) that is required so that the dam can fulfil all its purposes. Conversely, it can also
detect infeasibilities due to flood issues (i.e. too much water into the reservoir).

  * `solution`: a solution to evaluate for feasibility (implicitly, a reservoir). This solution is considered
    as a purpose
  * `inflow`: the inflow to test the solution against.
  * `forceInitialCondition`: imposes the value of the initial condition to be the worst case for the source
    of unfeasibility (i.e. start with a very low reservoir level if the infeasibility is due to water shortage,
    as full as possible otherwise).

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
function purposeShortage(solution::ReservoirSolution, inflow::Array{Float64, 1}, name::AbstractString="";
                         useOutputs::Array{Symbol, 1}=Symbol[], forceInitialCondition::Bool=false,
                         output="", verbosity::Int=1)
  # How should this model work? Lower the value of the purposes for each time step (say delta_t >= 0), minimise the sum
  # of all those deltas: if the objective value is zero, then there is a solution for this problem. Otherwise,
  # the objective value gives the shortage (from which indicators may be derived, such as a shortage index or
  # a vulnerability index).
  # Sub-zero objective values must be disallowed, otherwise the problem is likely to be unbounded.
  #     storage_t+1 = storage_t + inflow_t - (purpose_t - delta_t)
  #
  # A similar model may be used to detect infeasibility due to flood management (i.e. too much water flowing into the
  # reservoir): minimise the amount of water that is NOT in the inflow such that there is a solution.
  #     storage_t+1 = storage_t + (inflow_t - delta_t) - purpose_t
  #
  # The difference with the previous state equation is the sign of the delta_t variable; the objective is the same.
  # Nevertheless, both cannot be merged as one program: the solver will not choose "automagically" the worst case
  # between both approaches.
  # As a consequence, both can be merged as one program that allows delta_t to take any value (either positive or
  # negative). The absolute value of delta_t is then minimised, so that zero vaues are preferred over the other ones
  # and to avoid unboundedness.
  #     storage_t+1 = storage_t + inflow_t - purpose_t + delta_t

  if output != ""
    # Problem here: multiple models solved.
    error("Not implemented yet; needs a bit of refactoring everywhere before doing this properly")
  end

  if ! solution.feasible
    error("Received an infeasible solution: " * getVariant(getReservoir(solution)))
  end

  reservoir = getReservoir(solution)
  if verbosity >= 1
    println(" == ")
    println(" ==== Water shortage for purposes: ", getVariant(reservoir))
  end

  time_start = float(time_ns())

  ### Deal with input data.
  time_horizon = length(inflow)
  minLevel = toLength(getSolution(solution), time_horizon)
  period = getPeriod(reservoir)

  if verbosity >= 2
    println(" ====== Input time horizon: ", time_horizon)
  end

  ### Solve the identification optimisation models.
  # Not implemented as a maximum over the two options: rather, all tests are performed afterwards. This allows for
  # full flexibility.
  # First test for a shortage: if feasible for a zero objective, then the curve is feasible for the output; if feasible
  # for a nonzero objective, then reason found. If infeasible, then no conclusion.
  # Second test for a surplus: same kind of conclusions, but for a shortage.
  # Third test: if the the previous tests are not conclusive, then perform a third test that allows shortages and
  # surpluses at the same time.
  modelling_times = Float64[]
  writing_times = Float64[]
  solving_times = Float64[]

  (status_1, objective_1, delta_1, modelling_time_1, writing_time_1, solving_time_1) =
    __modelShortage(reservoir, time_horizon, forceInitialCondition, inflow, period, useOutputs, minLevel, verbosity)
  push!(modelling_times, modelling_time_1)
  push!(writing_times, writing_time_1)
  push!(solving_times, solving_time_1)

  if status_1 == :Optimal
    # Can get a conclusion from this first solve.
    if abs(objective_1) < 1.e-6 # No problem, this is feasible.
      delta = zeros(time_horizon)
      infeasibilityCause = :None
    else
      delta = delta_1
      infeasibilityCause = :Shortage
    end
  else
    # No conclusion for a shortage. Try a surplus.
    (status_2, objective_2, delta_2, modelling_time_2, writing_time_2, solving_time_2) =
      __modelSurplus(reservoir, time_horizon, forceInitialCondition, inflow, period, useOutputs, minLevel, verbosity)
    push!(modelling_times, modelling_time_2)
    push!(writing_times, writing_time_2)
    push!(solving_times, solving_time_2)

    if status_2 == :Optimal
      # Can get a conclusion from this second solve.
      if abs(objective_2) < 1.e-6 # No problem, this is feasible. (Should have been caught earlier, though.)
        delta = zeros(time_horizon)
        infeasibilityCause = :None
      else
        delta = delta_2
        infeasibilityCause = :Surplus
      end
    else
      # No conclusion for shortages or surpluses. Try a mix: always conclude.
      (status_3, objective_3, delta_3, modelling_time_3, writing_time_3, solving_time_3) =
        __modelMixed(reservoir, time_horizon, forceInitialCondition, inflow, period, useOutputs, minLevel, verbosity)
      push!(modelling_times, modelling_time_3)
      push!(writing_times, writing_time_3)
      push!(solving_times, solving_time_3)

      if status_3 == :Optimal
        # Determine whether it's mostly shortage or surplus.
        if abs(objective_3) < 1.e-6 # No problem, this is feasible. (Should have been caught earlier, though.)
          delta = zeros(time_horizon)
          infeasibilityCause = :None
        elseif objective_3 > 0 # Mostly shortage (if it was only shortage, it would have been caught earlier).
          delta = delta_3
          infeasibilityCause = :MostlyShortage
        else # Mostly surplus.
          delta = delta_3
          infeasibilityCause = :MostlySurplus
        end
      else
        error("Assertion exception: the mixed model could not give any answer.")
      end
    end
  end

  # Output if needed.
  if verbosity >= 1
    println(" ==== Infeasibility cause: ", string(infeasibilityCause), ".")
    println(" ==== Total shortage/surplus: ", sum(delta) * 1000000, " m^3 over ", time_horizon * getPeriod(reservoir), ".")
  end

  time_done = float(time_ns())

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
  options = ShortageOptions(useOutputs, forceInitialCondition)
  solution = ShortageSolution(infeasibilityCause, delta)
  time = ShortageSolverStatistics((time_done - time_start) / 1000000, modelling_times, writing_times, solving_times)
  solName = (name == "") ? Nullable{AbstractString}() : Nullable(name)
  return ReservoirSolution(reservoir, true, solution, solName, :ShortageSolver, "water shortage/surplus", options, time)
end
