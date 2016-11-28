["""
Evaluates a rule curve based on its hydropower potential.
"""]



immutable HydropowerOptions <: EvaluationOptions
  model::Symbol
  electricityPrices::Array{Float64, 1}
end

immutable HydropowerSolution <: EvaluationSolution
  totalRevenue::Float64
  revenues::Array{Float64, 1}
end

immutable HydropowerSolverStatistics <: EvaluationSolverStatistics
  total::Float64
end



"""
Estimate the hydropower potential of some solution with respect to a series of scenarios.

Inputs:

  * `solution`: a solution to evaluate for feasibility (implicitly, a reservoir)
  * `electricityPrices`: the eletricity prices to consider for evaluation (scenarios match those of the inflow).
    First indexed by scenario, then timestep, i.e. `electricityPrices[scenario, t]`
  * `model`: the model to use for evaluation (`:Nonconvex` or `:Convex`)
  * `probabilities`: probabilities for each scenario in the reservoir (default: `zeros(Float64, 0)`, indicating that the scenarios are equiprobable).
  * `imposeDiverted`: whether the inflow from the diverted rivers is imposed in the optimisation program, or left as a degree of freedom.
    This option is disabled by default: when the minimum for an objective function is reached, the diversions will be maximum
    (to minimise what must be stored).

To evaluate a given solution, the search space can be limited by rule curves:

  * `minLevel`: impose a minimum rule curve to follow (usually used with `inflow` to evaluate a solution)
  * `maxLevel`: impose a maximum rule curve to follow (usually used with `inflow` to evaluate a solution)

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

TODO: partially PWL model? (Height-volume and/or objective function.) See PiecewiseLinearReservoirBathymetry (to complete).
TODO: What about price uncertainty? Should be enough to just take the lowest ones… (robust approach).
TODO: Power depends on scenario? Should be configurable! If it actually depends: evaluate the POTENTIAL for hydropower, as the best decisions are always taken. Otherwise: how are dams managed for hydropower? May need to go for multistage stochastic things if planned one month in advance, or nothing can be done if it is market-controlled.
TODO: only one hydropower unit considered in the model.
TODO: minimum head to produce power? Allowable to go below (but lose production)?
"""
function hydropower(solution::ReservoirSolution, electricityPrices::Array{Float64, 1}, model::Symbol=:nonconvex;
                    probabilities::Array{Float64,1}=zeros(Float64, 0), imposeDiverted::Bool=true, inflow::Array{Float64, 2}=zeros(Float64, 0, 0),
                    minLevel=EmptyInverseSolution(), maxLevel::InverseSolution=EmptyInverseSolution(),
                    useOutputs::Array{Symbol, 1}=Symbol[],
                    name="", verbosity::Int=1)
  reservoir = getReservoir(solution)

  # TODO: if hydropower not in useOutputs, shout, scream, yell.

  # TODO: this is awkwardly strange, think again about those minLevel/maxLevel/solution to evaluate. Actually, here, the solution to evaluate is the minimum curve--for now.
  if minLevel == EmptyInverseSolution()
    minLevel = solution
  else
    error("Currently, the code does strange things. Please excuse.")
  end

  recognised_models = [:Nonconvex, :Convex]
  if ! in(model, recognised_models)
    error("Unrecognised model: " * string(model) * ". Please choose from the list: " * string(recognised_models))
  end

  if isempty(bathymetry(reservoir)) || isblackbox(bathymetry(reservoir))
    error("Hydropower evaluation requires full reservoir bathymetry information!")
  end

  if length(hydropower(reservoir)) < 1
    error("Hydropower evaluation requires at least one hydropower unit!")
  end

  if length(hydropower(reservoir)) != 1
    error("TODO: only one hydropower unit supported for now!")
  end

  if verbosity >= 1
    println(" == ")
    println(" ==== Hydropower evaluation for the reservoir " * getName(reservoir) * ", scenario generation " * getVariant(reservoir))
  end

  ## Define some constants.
  ρ = 1000 # Water density [kg/m^3]
  g = 9.81 # Earth gravity [m/s²]

  time_start = float(time_ns())

  ## Deal with input data.
  # TODO: almost copy-paste from rulecurve_min_default.
  forced_inflow = ! isempty(inflow)
  forced_min = ! isempty(minLevel)
  forced_max = ! isempty(maxLevel)

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

  ## Build the optimisation models.
  # Two main variants: either the direct, nonconvex formulation (through JuMP), or a convexified version thereof (through Convex).
  if model == :Nonconvex
    # Build a JuMP, nonconvex model.
    # Based on doi:10.1061/(ASCE)0733-9496(2003)129%3A3(178), with a single objective function (revenue).
    m = Model(solver=GurobiSolver(OutputFlag=(verbosity >= 4)))

    @variable(m, getMinimumCapacity(reservoir) <= storage[1:scenarios, 0:time_horizon] <= getMaximumCapacity(reservoir)) # 10^6 m^3
    @variable(m, release[1:scenarios, 1:time_horizon] >= 0) # 10^6 m^3/time step
    @variable(m, 0 <= releaseHydropower[1:scenarios, 1:time_horizon] <= getMaximumDischarge(getHydropower(reservoir))) # 10^6 m^3/time step
    @variable(m, 0 <= hydropower[1:scenarios, 1:time_horizon] <= getMaximumPower(getHydropower(reservoir)) / (ρ * g * efficiency(re))) # W / ((kg / m^3) * (m / s²) * 1) (power considered constant over the time step, so equivalent to energy for optimisation)
    @variable(m, minCapacity(reservoir) <= storage[1:scenarios, 0:time_horizon] <= maxCapacity(reservoir)) # 10^6 m^3
    @variable(m, level(reservoir, minCapacity(reservoir)) <= level[1:scenarios, 0:time_horizon] <= computeLevel(reservoir, maxCapacity(reservoir))) # 10^6 m^3
    if ! imposeDiverted
      @variable(m, divert[1:scenarios, 1:time_horizon, diverted_names] >= 0) # 10^6 m^3/time step.
    end

    @objective(m, Min, sum{probabilities[scenario] * electricityPrices[scenario, t] * hydropower[scenario, t], scenario in 1:scenarios, t in 1:time_horizon})

    # Define the hydropower variable.
    @NLconstraint(m, c_hydropower[scenario in 1:scenarios, t in 1:time_horizon], hydropower[scenario, t] == level[scenario, t] * releaseHydropower[scenario, t])

    for scenario = 1:scenarios
      for t = 1:time_horizon
        # Define the reservoir level based on the volume.
        hlp_JuMP_bathymetry(m, storage, level, reservoir, nonlinear_mode=:Nonconvex)

        # Deal with the inflow: either imposed via a parameter, or computed based on the rivers information.
        if ! forced_inflow
          natural_inflow = sum([getScenario(r, scenario, t) for r in naturalRivers(reservoir)])

          # Compute the potential diverted inflow to the reservoir.
          if imposeDiverted
            diverted_inflow = sum([maxAllowableFlow(r, scenario, t) for r in divertedRivers(reservoir)])
          else
            diverted_inflow = sum([divert[scenario, t, rn] for rn in diverted_names])
            for r in divertedRivers(reservoir)
              @constraint(m, divert[scenario, t, getName(r)] <= maxAllowableFlow(r, scenario, t))
            end
          end

          total_inflow = natural_inflow + diverted_inflow
        else
          total_inflow = inflow[scenario, t]
        end

        # Implement the state equation.
        @constraint(m, storage[scenario, t] == storage[scenario, t - 1] - totalPurposes(reservoir, period) + total_inflow - release[scenario, t] - release[scenario, t])
        if ! isinf(max_releases); @constraint(m, release[scenario, t] <= max_releases); end
      end
    end

    # Implement the rule curves if given.
    if forced_min
      @constraint(m, storage[scenario, t] >= minLevel[t])
    end
    if forced_max
      @constraint(m, storage[scenario, t] <= maxLevel[t])
    end
  elseif model == :Convex
    # Build a Convex model.
    error("TODO")
  else
    error("Non recognised model: " * string(model))
  end

  if length(probabilities) == 0
    probabilities = ones(Float64, scenarios) / scenarios
    if verbosity >= 2
      println(" ====== Defining probabilities for scenarios: equiprobable")
    end
  end

  # Adapt the array sizes of the maximum/minimum rule curves, if needed.
  if forced_min
    minLevel = getSolution(minLevel)
    if length(minLevel) > time_horizon
      minLevel = minLevel[1:time_horizon]
    elseif length(minLevel) < time_horizon
      minLevel = repmat(minLevel, 1, Int(round(time_horizon / length(minLevel))))
    end
  end
  if forced_max
    maxLevel = getSolution(maxLevel)
    if length(maxLevel) > time_horizon
      maxLevel = maxLevel[1:time_horizon]
    elseif length(maxLevel) < time_horizon
      maxLevel = repmat(maxLevel, 1, Int(round(time_horizon / length(maxLevel))))
    end
  end

  # Complete the process.
  time_end = float(time_ns())
  total_proba /= folds

  # A bit of printing.
  if verbosity >= 1
    println(" ==== Probability of being infeasible: ", total_proba * 100, "%.")
    println(" ==== Total time: ", (time_end - time_start) / 1000000, "ms")
    if verbosity >= 2
      println(" ====== Time to generate LS and TS: ", total_time_split / 1000000, "ms")
      println(" ====== Time to optimise: ", total_time_solved / 1000000, "ms")
      println(" ====== Time to check feasibility: ", total_time_evaluated / 1000000, "ms")
    end
    println(" == ")
  end

  # Prepare the output.
  options = HydropowerOptions(solver, sizes, folds)
  solution = HydropowerSolution(total_proba)
  time = HydropowerSolverStatistics((time_end - time_start) / 1000000, total_time_split / 1000000, total_time_solved / 1000000, total_time_evaluated / 1000000)
  solName = (name == "") ? Nullable{AbstractString}() : Nullable(name)
  return ReservoirSolution(reservoir, true, solution, solName, :HydropowerSolver, "Cross-validation feasibility estimation", options, time) # TODO: feasibility "true" does not make sense here!
end
