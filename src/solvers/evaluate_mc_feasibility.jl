["""
Evaluates a rule curve based on a Monte Carlo approach to estimate its feasibility probability.
"""]



immutable FeasibilityMonteCarloOptions <: EvaluationOptions
  performScenarioGeneration::Bool
  method::Symbol
  numberScenarios::Int # Number of scenarios, either generated or given to the evaluator.
end

immutable FeasibilityMonteCarloSolution <: EvaluationSolution
  numberScenarios::Int
  countInfeasible::Int
  countInfeasibleSurplus::Int
  countInfeasibleShortage::Int

  probabilityInfeasible::Float64
  probabilityInfeasibleSurplus::Float64
  probabilityInfeasibleShortage::Float64

  shortage::Array{Float64, 2}
  shortageIndex::Array{Float64, 1}
  shortageVulnerability::Array{Float64, 1}

  averageShortage::Float64
  averageShortageIndex::Float64
  averageShortageVulnerability::Float64
end

immutable FeasibilityMonteCarloSolverStatistics <: EvaluationSolverStatistics
  total::Float64
  generation::Float64
  check::Float64
end



"""
Estimate the feasibility probabilities of the given solution with a Monte-Carlo approach, generating
many scenarios and testing the solution against each of them. Testing the solution means that the solver
detects whether there is some infeasibility, and in that case, determines its cause (surplus or shortage of water).

Inputs:

  * `solution`: a solution to evaluate for feasibility (implicitly, a reservoir)
  * `method`: the scenario generating method to use (see statistical module, function `generateScenarios`);
    this parameter is ignored if `scenarios` is set.
  * `numberScenarios`: the number of scenarios to generate in order to perform the estimation;
    this parameter is ignored if `scenarios` is set.
  * `scenarioDuration`: the duration for each generated scenario; this parameter is ignored if `scenarios` is set.
  * `scenarios`: allows to bypass the scenario generation step; `eachscenario(scenarios)` must iterate
    through the requested scenarios `s` and `totalInflow`. Scenario generation is still performed if
    `isempty(scenarios)` is `true`.
    The parameters `method`, `numberScenarios`, and `scenarioDuration` are all ignored if this is set.

Finally, `verbosity` controls the level of output from this solver:

  * `-1`: no output
  * `0`: no output, except for solver's warnings
  * `1`: summary output
  * `2`: detailed output of this algorithm's details
  * `3`: detailed output of this algorithm's and its sub-algorithms' details
  * `4`: including solver output
"""
function feasibilityProbability(solution::ReservoirSolution;
                                method::Symbol=:Adam2015, numberScenarios::Int=1000, scenarioDuration::Period=Year(1),
                                scenarios=zeros(Float64, 0), useOutputs::Array{Symbol, 1}=Symbol[],
                                name="", verbosity::Int=1)
  if verbosity >= 1
    println(" == ")
    println(" ==== Monte-Carlo simulation to estimate feasibility of solution: ", getVariant(solution))
    if ! isempty(scenarios)
      if applicable(getVariant, scenarios)
        println(" ==== Testing against the given set of scenarios: ", getVariant(scenarios))
      else
        println(" ==== Testing against the given set of scenarios. ")
      end
    end

    if verbosity >= 2
      if isempty(scenarios)
        println(" ====== Going to perform " * string(numberScenarios) * " simulations, each of them for a duration of " * string(scenarioDuration))
        println(" ====== Method to generate scenarios: " * string(method))
      else
        println(" ====== Bypassed scenario generation")
        println(" ====== Going to perform " * string(length(eachscenario(scenarios))) * " simulations")
      end
    end
  end

  reservoir = solution.reservoir
  time_start = float(time_ns())

  ## Generate scenarios if need be.
  mustGenerateScenarios = isempty(scenarios)
  if mustGenerateScenarios
    scenarios = generateScenarios(reservoir, numberScenarios, getPeriod(reservoir), method, scenarioDuration, iterator=true)
  end
  n_scenarios = countScenarios(scenarios)
  n_timesteps = countTimeSteps(scenarios)
  time_gen = float(time_ns())

  if n_scenarios < 1
    error("Zero scenarios given, but scenario generation not disabled. Please ensure that isempty(scenarios) when calling the function with no scenarios.")
  end

  ## Actual probability computation.
  count = 0
  countInfeasibleSurplus = 0
  countInfeasibleShortage = 0
  shortages = zeros(Float64, n_scenarios, n_timesteps)
  shortageIndices = zeros(Float64, n_scenarios)
  shortageVulnerabilities = zeros(Float64, n_scenarios)
  for s in eachscenario(scenarios)
    count += 1
    total_s = totalInflow(s)
    total_s = reshape(total_s, 1, length(total_s))

    sol = purposeShortage(solution, vec(total_s), useOutputs=useOutputs, forceInitialCondition=false, verbosity=(verbosity >= 3) ? verbosity : -1)
    if ! isFeasible(sol)
      println(" ====== /!\\ Evaluation model failed! ")
    else
      if in(getSolutionObject(sol).infeasibilityCause, [:Shortage, :MostlyShortage])
        countInfeasibleShortage += 1
        shortages[count, :] = getSolutionObject(sol).delta # Count only actual shortage here.
      elseif in(getSolutionObject(sol).infeasibilityCause, [:Surplus, :MostlySurplus])
        countInfeasibleSurplus += 1
      end
      shortageVulnerabilities[count] = computePurposeVulnerabilityIndex(sol, getPurposes(reservoir))
      shortageIndices[count] = computePurposeShortageIndex(sol, getPurposes(reservoir))
    end
  end
  time_check = float(time_ns())

  ## Prepare output.
  countInfeasible = countInfeasibleShortage + countInfeasibleSurplus
  probabilityInfeasible = countInfeasible / count
  meanShortages = mean(shortages)
  meanShortageIndices = mean(shortageIndices)
  meanShortageVulnerabilities = mean(shortageVulnerabilities)

  if verbosity >= 1
    println(" ==== Probability of being infeasible: ", probabilityInfeasible * 100, "%.")
    if verbosity >= 2
      println(" ====== Total number of trials: " * string(count))
      println(" ====== Total number of infeasibilities: " * string(countInfeasible))
      println(" ====== Total number of infeasibilities due to water shortage: " * string(countInfeasibleShortage))
      println(" ====== Total number of infeasibilities due to water surplus: " * string(countInfeasibleSurplus))
      println(" ====== Average shortage (10^6 m^3): " * string(meanShortages))
      println(" ====== Average shortage index: " * string(meanShortageIndices))
      println(" ====== Average shortage vulnerability: " * string(meanShortageVulnerabilities))
    end
  end

  if verbosity >= 1
    println(" ==== Total time: ", (time_check - time_start) / 1000000, "ms")
    if verbosity >= 2
      mustGenerateScenarios && println(" ====== Time to generate the scenarios (creating an iterator): ", (time_gen - time_start) / 1000000, "ms")
      println(" ====== Time to check individual scenarios: ", (time_check - time_gen) / 1000000, "ms")
    end
    println(" == ")
  end

  options = FeasibilityMonteCarloOptions(mustGenerateScenarios, method, count)
  solution = FeasibilityMonteCarloSolution(count, countInfeasible, countInfeasibleSurplus, countInfeasibleShortage,
                                           probabilityInfeasible, countInfeasibleSurplus / count, countInfeasibleShortage / count,
                                           shortages, shortageIndices, shortageVulnerabilities,
                                           meanShortages, meanShortageIndices, meanShortageVulnerabilities)
  time = FeasibilityMonteCarloSolverStatistics((time_check - time_start) / 1000000, (time_gen - time_start) / 1000000, (time_check - time_gen) / 1000000)
  solName = (name == "") ? Nullable{AbstractString}() : Nullable(name)
  return ReservoirSolution(reservoir, true, solution, solName, :FeasibilityMonteCarloSolver, "Monte Carlo feasibility estimation", options, time)
end
