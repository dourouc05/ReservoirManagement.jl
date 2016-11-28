["""
The most basic solver: optimises a minimum rule curve over scenarios.
"""]



immutable DefaultMinimumRuleCurveOptions <: MinimumRuleCurveOptions
  scenarioProbabilities::Array{Float64, 1}
  objective::Symbol
  imposeDiverted::Bool
  useOutputs::Array{Symbol, 1}
end

"""
The default solver solution has a rule curve and the solution to each of the scenarios that are given in input.
"""
immutable DefaultMinimumRuleCurveSolution <: MinimumRuleCurveSolution
  solution::Array{Float64, 1}
  levels::Array{Float64, 2}
end

countTimeSteps(sol::DefaultMinimumRuleCurveSolution) = length(sol.solution)

immutable DefaultMinimumRuleCurveSolverStatistics <: MinimumRuleCurveSolverStatistics
  total::Float64
  modelling::Float64
  writingModelFile::Float64
  solving::Float64
  postprocessing::Float64
end



"""
Starts a configurable optimisation program on the given water reservoir.

  * `reservoir`: the reservoir to compute a rule for.
  * `probabilities`: probabilities for each scenario in the reservoir
    (default: `zeros(Float64, 0)`, indicating that the scenarios are equiprobable).
  * `objective`: :Total (minimise the total storage of the reservoirs among all scenarios)
    or :Rule (minimise the curve above all scenarios such that feasibility is guaranteed).
  * `imposeDiverted`: whether the inflow from the diverted rivers is imposed in the optimisation program,
    or left as a degree of freedom.
    This option is disabled by default: when the minimum for an objective function is reached, the diversions
    will be maximum(to minimise what must be stored).
  * `inflow`: impose exterior inflow for the model (parts of `reservoir` are then ignored)

To evaluate a given solution, the search space can be limited by rule curves; these bounds are usually used with
the `inflow` parameter.

  * `minLevel`: impose a minimum rule curve to follow
  * `maxLevel`: impose a maximum rule curve to follow

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
function minimumRuleCurve(reservoir::Reservoir, name::AbstractString="";
                          probabilities::Array{Float64,1}=zeros(Float64, 0), objective::Symbol=:Rule, imposeDiverted::Bool=true, inflow::Array{Float64, 2}=zeros(Float64, 0, 0),
                          minLevel=EmptyInverseSolution(), maxLevel::InverseSolution=EmptyInverseSolution(),
                          useOutputs::Array{Symbol, 1}=Symbol[], outputLength::Int=countTimeSteps(reservoir),
                          output="", verbosity::Int=1)
  recognised_objectives = [:Total, :Rule]
  if ! in(objective, recognised_objectives)
    error("Unrecognised objective: " * string(objective) * ". Please choose from the list: " * string(recognised_objectives))
  end

  if verbosity >= 1
    println(" == ")
    println(" ==== Optimisation model for operational rules: ", getVariant(reservoir))
  end

  time_start = float(time_ns())

  ## Deal with input data.
  forced_inflow = ! isempty(inflow)
  forced_min = ! isempty(minLevel)
  forced_max = ! isempty(maxLevel)

  period = getPeriod(reservoir)
  scenarios = forced_inflow ? size(inflow, 1) : countScenarios(reservoir)
  time_horizon = forced_inflow ? size(inflow, 2) : countTimeSteps(reservoir)
  diverted_names = getDivertedRiverNames(reservoir)

  if verbosity >= 2
    if forced_inflow
      println(" ====== Inflow takes forced values (not taken from input reservoir)")
    end
    if forced_min
      println(" ====== Minimum level imposed from a minimum rule curve")
    end
    if forced_max
      println(" ====== Maximum level imposed from a maximum rule curve")
    end
    println(" ====== Input number of scenarios: ", scenarios)
    println(" ====== Input time horizon: ", time_horizon)
  end

  if length(probabilities) == 0
    probabilities = ones(Float64, scenarios) / scenarios
    if verbosity >= 2
      println(" ====== Defining probabilities for scenarios: equiprobable")
    end
  end

  # Adapt the array sizes of the maximum/minimum rule curves, if needed.
  if forced_min
    minLevel = toLength(getSolution(minLevel), time_horizon)
  end
  if forced_max
    maxLevel = toLength(getSolution(maxLevel), time_horizon)
  end

  ## Create the optimisation model.
  # Two main variants of the model: let the solver decide for the diversions (imposeDiverted == false), or impose it to the maximum (imposeDiverted == true).
  # The solver should use this degree of freedom as much as possible to lower the need for storage, i.e. use this diverted river a lot.
  m = Model(solver=GurobiSolver(OutputFlag=(verbosity >= 4)))
  @variable(m, minCapacity(reservoir) <= storage[1:scenarios, 0:time_horizon] <= maxCapacity(reservoir)) # 10^6 m^3
  @variable(m, release[1:scenarios, 1:time_horizon] >= 0) # 10^6 m^3/time step
  if objective == :Rule
    @variable(m, minCapacity(reservoir) <= rule[1:time_horizon] <= maxCapacity(reservoir)) # 10^6 m^3
  end
  if ! imposeDiverted
    @variable(m, divert[1:scenarios, 1:time_horizon, diverted_names] >= 0) # 10^6 m^3/time step.
  end

  # Objective function.
  if objective == :Rule
    @objective(m, Min, sum(rule))
  elseif objective == :Total
    @objective(m, Min, sum{ probabilities[scenario] * sum(storage[scenario, :]), scenario = 1:scenarios })
  end

  # Main set of constraints.
  for scenario = 1:scenarios
    for t = 1:time_horizon
      # Deal with the inflow: either imposed via a parameter, or computed based on the rivers information.
      if ! forced_inflow
        if imposeDiverted
          total_inflow = hlp_JuMP_inflow_discharge(m, reservoir, scenario, t)
        else
          total_inflow = hlp_JuMP_inflow_discharge(m, reservoir, scenario, t; imposeDiverted=divert)
        end
      else
        total_inflow = inflow[scenario, t]
      end

      # Implement the state equation.
      @constraint(m, storage[scenario, t] == storage[scenario, t - 1] - release[scenario, t]
                    - hlp_JuMP_purposes_discharge(m, reservoir, period) + total_inflow)
      hlp_JuMP_outputs_constraint(m, release[scenario, t], storage[scenario, t], reservoir, useOutputs)

      # Implement the rule curves if given.
      if forced_min
        @constraint(m, storage[scenario, t] >= minLevel[t])
      end
      if forced_max
        @constraint(m, storage[scenario, t] <= maxLevel[t])
      end
    end
  end

  # Other constraints, depending on the parameters.
  if objective == :Rule
    for scenario = 1:scenarios
      for t = 1:time_horizon
        @constraint(m, storage[scenario, t] <= rule[t])
      end
    end
  end

  time_modelled = float(time_ns())

  ## Services around the model.
  if output != ""
    writeLP(m, output)
  end
  time_output = float(time_ns())

  status = solve(m, suppress_warnings=verbosity < 0); time_solved = float(time_ns())
  if status != :Infeasible
    # levels = getvalue(storage[:, 1:end])
    # The following is a workaround for JuMP 0.12+: https://github.com/JuliaOpt/JuMP.jl/issues/730 (keyword end does not work).
    levels = getvalue(storage)[:, 1:time_horizon]

    if objective == :Rule
      target = getvalue(rule)
    elseif objective == :Total
      target = vec(maximum(levels, 1))
    end

    if verbosity >= 1
      println(" ==== Total storage: ", getobjectivevalue(m) * 1000000, " m^3 over ", time_horizon * period, ".")
    end

    target = target[1:end] # Remove time step 0 (i.e. initial condition).
    feasible = true
  else
    if verbosity >= 1
      println(" ==== Infeasible!")
    end
    feasible = false
    target = zeros(time_horizon) / 0. # NaN as 0/0.
    levels = zeros(time_horizon, scenarios) / 0. # NaN as 0/0.
  end

  time_done = float(time_ns())

  # Compute the timings (in ms).
  total_time = (time_done - time_start) / 1000000
  modelling_time = (time_modelled - time_start) / 1000000
  writing_time = (time_output - time_modelled) / 1000000
  solving_time = (time_solved - time_output) / 1000000
  postprocessing_time = (time_done - time_solved) / 1000000

  # Final output.
  if verbosity >= 1
    println(" ==== Total time: ", total_time, "ms")
    if verbosity >= 2
      println(" ====== Time to create the JuMP model: ", modelling_time, "ms")
      println(" ====== Time to write the LP model if required: ", writing_time, "ms")
      println(" ====== Time to actually solve the model: ", solving_time, "ms")
      println(" ====== Time to prepare the return values: ", postprocessing_time, "ms")
    end
    println(" == ")
  end

  # Prepare the data structures for the outside world.
  options = DefaultMinimumRuleCurveOptions(probabilities, objective, imposeDiverted, useOutputs)
  solution = DefaultMinimumRuleCurveSolution(target[1:outputLength], levels)
  time = DefaultMinimumRuleCurveSolverStatistics(total_time, modelling_time, writing_time, solving_time, postprocessing_time)
  solName = (name == "") ? Nullable{AbstractString}() : Nullable(name)
  return ReservoirSolution(reservoir, feasible, solution, solName, :BasicRuleCurveSolver, "basic rule curve from inflow", options, time)
end
