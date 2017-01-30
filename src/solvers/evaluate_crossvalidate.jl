["""
Evaluates a rule curve based on a cross-validation approach to esimate feasibility probabilities.
"""]



immutable FeasibilityCrossValidationOptions <: EvaluationOptions
  solver::Function
  sizes::Tuple{Float64, Float64}
  folds::Int
end

immutable FeasibilityCrossValidationSolution <: EvaluationSolution
  numberFolds::Int

  probabilityInfeasible::Float64
  probabilityInfeasibleSurplus::Float64
  probabilityInfeasibleShortage::Float64
  probabilityNoRuleCurve::Float64

  averageShortage::Float64
  averageShortageIndex::Float64
  averageShortageVulnerability::Float64
end

immutable FeasibilityCrossValidationSolverStatistics <: EvaluationSolverStatistics
  total::Float64
  split::Float64
  optimisation::Float64
  check::Float64
end



"""
Estimate the feasibility probabilities of the given solution with cross-validation. This function handles
the whole process, from a separation in a learning set (LS) and a test set (TS), to computing a solution
with the learning set, to validate it against the test set. The division into LS and TS is performed randomly
for each fold using random subsampling.

Inputs:

  * `reservoir`: a reservoir to evaluate the given method on
  * `function`: the solver that computes the solution to evaluate. Its only inputs must be the reservoir
    to optimise (type: `Reservoir`) and the verbosity, while the output must be an untouched `ReservoirSolution`.
  * `sizes`: the relative sizes of the learning set and the test set. The two values must sum up to 1.
  * `folds`: the number of times the total process should be repeated; the output values are then
    averaged over all folds (repetitions).
  * `verbosity`: level of output.
    * `-1`: no output
    * `0`: no output, except for solver's warnings
    * `1`: summary output
    * `2`: detailed output of this algorithm's details
    * `3`: detailed output of this algorithm's and its sub-algorithms' details
    * `4`: including solver output

TODO: what about solver identification in the output data structure?
"""
function crossValidatedFeasibilityProbability(reservoir::Reservoir, solver::Function;
                                              sizes::Tuple{Float64, Float64}=(.8, .2),
                                              folds::Int=1, guaranteePeriods::Int=countTimeSteps(reservoir),
                                              name="", verbosity::Int=1)
  @assert folds >= 1 "At least one repetition"
  @assert countScenarios(reservoir) >= 2 "At least two scenarios are required in order to make two sets."

  if verbosity >= 1
    println(" == ")
    println(" ==== Cross-validation procedure for the reservoir " * getName(reservoir) * ", scenario generation " * getVariant(reservoir))
    if verbosity >= 2
      println(" ====== Going to perform " * string(folds) * " repetition" * ((folds > 1) ? "s" : ""))
    end
  end

  time_start = float(time_ns())
  total_time_split = 0.
  total_time_solved = 0.
  total_time_evaluated = 0.
  total_proba_inf = 0.
  total_proba_sur = 0.
  total_proba_srt = 0.
  total_no_rule = 0.
  total_shortage = 0.
  total_shrt_idx = 0.
  total_shrt_vul = 0.

  for k = 1:folds
    if verbosity >= 1 && k % 10 == 0; println(" ====== Current fold: " * string(k)); end

    time_iter = float(time_ns())

    # Generate the various sets. Limit the test set to the guarantee length
    if verbosity >= 2
      println(" ====== Splitting the scenarios in a learning set and a test set. ")
    end
    split = splitIntoSets(reservoir, collect(sizes))
    LS = split[1]
    TS = keepTimeSteps(split[2], 1:guaranteePeriods)
    time_split = float(time_ns())

    # Learn on the learning set.
    if verbosity >= 2
      println(" ====== Optimising a solution on the learning set. ")
    end
    solution = solver(LS, verbosity >= 3 ? verbosity : -1)
    time_solved = float(time_ns())

    if solution.feasible
      # Test on the test set.
      if verbosity >= 2
        println(" ====== Evaluating the solution on the test set. ")
      end
      evaluation = feasibilityProbability(solution, scenarios=TS, verbosity=verbosity >= 3 ? verbosity : -1)
      time_evaluated = float(time_ns())

      # Update the global counters.
      total_time_split += time_split - time_iter
      total_time_solved += time_solved - time_split
      total_time_evaluated += time_evaluated - time_solved
      total_proba_inf += evaluation.solution.probabilityInfeasible
      total_proba_sur += evaluation.solution.probabilityInfeasibleSurplus
      total_proba_srt += evaluation.solution.probabilityInfeasibleShortage
      # total_no_rule
      total_shortage += evaluation.solution.averageShortage
      total_shrt_idx += evaluation.solution.averageShortageIndex
      total_shrt_vul += evaluation.solution.averageShortageVulnerability
    else
      # No rule curve could be derived. Update the global counters.
      total_time_split += time_split - time_iter
      total_time_solved += time_solved - time_split
      # total_time_evaluated
      # total_proba_inf
      # total_proba_sur
      # total_proba_srt
      total_no_rule += 1
      # total_shortage
      # total_shrt_idx
      # total_shrt_vul
    end
  end

  # Complete the process.
  time_end = float(time_ns())
  total_proba_inf /= folds
  total_proba_sur /= folds
  total_proba_srt /= folds
  total_no_rule /= folds
  total_shortage /= folds
  total_shrt_idx /= folds
  total_shrt_vul /= folds

  # A bit of printing.
  if verbosity >= 1
    println(" ==== Probability of being infeasible: ", total_proba_inf * 100, "%.")
    println(" ==== Probability of being infeasible due to water surplus:  ", total_proba_sur * 100, "%.")
    println(" ==== Probability of being infeasible due to water shortage: ", total_proba_srt * 100, "%.")
    println(" ==== Probability of having no rule curve: ", total_no_rule * 100, "%.")
    println(" ==== Average shortage (10^6 m^3): ", total_shortage)
    println(" ==== Average shortage index: ", total_shrt_idx)
    println(" ==== Average shortage vulnerability: ", total_shrt_vul)
    println(" ==== Total time: ", (time_end - time_start) / 1000000, "ms")
    if verbosity >= 2
      println(" ====== Time to generate LS and TS: ", total_time_split / 1000000, "ms")
      println(" ====== Time to optimise: ", total_time_solved / 1000000, "ms")
      println(" ====== Time to check feasibility: ", total_time_evaluated / 1000000, "ms")
    end
    println(" == ")
  end

  # Prepare the output.
  options = FeasibilityCrossValidationOptions(solver, sizes, folds)
  solution = FeasibilityCrossValidationSolution(folds, total_proba_inf, total_proba_sur, total_proba_srt, total_no_rule, total_shortage, total_shrt_idx, total_shrt_vul)
  time = FeasibilityCrossValidationSolverStatistics((time_end - time_start) / 1000000, total_time_split / 1000000, total_time_solved / 1000000, total_time_evaluated / 1000000)
  solName = (name == "") ? Nullable{AbstractString}() : Nullable(name)
  return ReservoirSolution(reservoir, true, solution, solName, :FeasibilityCrossValidationSolver, "Cross-validation feasibility estimation", options, time) # TODO: feasibility "true" does not make sense here!
end
