abstract MonteCarloScenarios
immutable FixedMonteCarloScenarios <: MonteCarloScenarios
  scenarios::Reservoir
end
immutable GeneratedMonteCarloScenarios <: MonteCarloScenarios
  method::Symbol
end

abstract ReportEvaluator
immutable CrossValidatedFeasibilityProbability <: ReportEvaluator end
immutable InverseModel <: ReportEvaluator end
immutable MonteCarloFeasibilityProbability <: ReportEvaluator
  scenarios::MonteCarloScenarios
end

"""
Prints a complete solution report for a reservoir, with multiple solutions and their evaluation.

Mandatory arguments:

  * `reservoir`: the reservoir to study
  * `ruleDuration`: the duration of the rule curves to produce
  * `guarantee`: the duration to guarantee for those rule curves

Optional arguments for input data:

  * `currentRuleCurve`: the existing rule curve, used as a benchmark for the newly computed ones.
    Its periodicity must be that of `reservoir`.
  * `climateChangeScenarios`: a set of inflow scenarios for some climate change, associating a description
    to the inflow. There is no distinction between the rivers in this argument!

Optional arguments to select the solvers and evaluators to use:

  * `scenarioGenerationAlgorithms`: algorithms used to generate scenarios (see Scenario generation module, the symbols
    are the argument to `scenarioGeneration`).
    Default: `[:Merge, :Mix]`
  * `ciAlgorithmsLevels`: algorithms used to compute the confidence intervals, associated with their confidence values
    (see Statistical module, the symbols are the argument to `generateScenarios`).
    Default: `Dict(:TStudent => [.95])`
  * `solvers`: solvers that are used to generate solutions (see Solver module, the symbols correspond to the function
    names of the solvers).
    Default: `[:minimumRuleCurve, :safeMinimumRuleCurve]`
  * `evaluators`: algorithms that are used to evaluate the new rule curves (see Solver module, the symbols correspond
    to the function names of the solvers).
    Default: `[]`
    TODO: not really flexible (think of parameters).

    The current list of possible evaluators is:

    * `:feasibilityProbabilitySyntheticScenarios`: computes the feasibility probability on synthetic scenarios
      (generated based on the statistical analysis of the data sets)
    * `:feasibilityProbabilityAllDatasets`: computes the feasibility probability on historical scenarios (all the given
       data sets), i.e. resubstitution
    * `:feasibilityProbabilityClimateChange`: computes the feasibility probability on provided climate change scenarios
      (where all rivers are merged)
    * `:hydropower`: evaluates the hydropower potential (with a generated price scenario)
      TODO: price scenario from outside.
    * `:worstCaseFunctions`:
    * `:maximumDeviation`:

  * `globalEvaluators`: algorithms that are used to evaluate algorithms that generate the rule curves, and not the
    rule curves themselves (see Solver module, the symbols correspond to the function names of the solvers).
    Default: `[]`
    TODO: not really flexible (think of parameters).
  * `disabledCombinations`: combinations of the above parameters that should *not* be run. This list contains pairs
    of a scenario generation algorithm and a solver. It is not possible to disable confidence intervals this way.
    Default: `[]`.
    Recommended value: `[(:Mix, :safeMinimumRuleCurve)]`, as this combination takes a lot of time to compute.
  * `verbosity`: passed as is to the solvers (`-1` to disable all output)
  * `saveFolder`: folder name to use with the following options
  * `savePlots`: whether the script should prepare some plots
  * `saveIntermediate`: whether the script should save the outputs for future reuse

TODO: how to add new solvers/evaluators? Work with a hash table and store lambdas with a fixed signature?
TODO: how to deal with solvers/evaluators parameters? Still via lambdas? Within the hash table?
TODO: provide more types to avoid typing Tuple{Tuple} too often. Use the defined data structures.
"""
function report(reservoir::Reservoir;
                ruleDuration::Period=Year(1), guarantee::Period=Year(2), ciPeriod::Period=Week(1),
                currentRuleCurve::Array{Float64, 1}=zeros(Float64, 0),
                climateChangeScenarios::Dict{Symbol, Array{TimeSeries.TimeArray{Float64,1,DateTime,Array{Float64,1}}, 1}}=Dict{Symbol, Array{TimeSeries.TimeArray{Float64,1,DateTime,Array{Float64,1}}, 1}}(),

                scenarioGenerationAlgorithms::Array{Symbol, 1}=Symbol[:Merge, :Mix],
                ciAlgorithmsLevels::Dict{Symbol, Array{Float64, 1}}=Dict{Symbol, Array{Float64, 1}}(:TStudent => [.95, .99]),
                solvers::Array{Symbol, 1}=[:minimumRuleCurve, :safeMinimumRuleCurve],
                evaluators::Array{Symbol, 1}=Symbol[],#:feasibilityProbabilityAllDatasets,:feasibilityProbabilityClimateChange], # :feasibilityProbability, # , :hydropower # , :worstCaseFunctions, :maximumDeviation
                globalEvaluators::Array{Symbol, 1}=Symbol[],#:crossValidatedFeasibilityProbability],
                disabledCombinations::Array{Tuple{Symbol, Symbol}, 1}=Tuple{Symbol, Symbol}[],#(:Mix, :safeMinimumRuleCurve)],
                saveFolder::AbstractString="", savePlots::Bool=false, saveIntermediate::Bool=false,
                verbosity::Int=1)
  # Some assumptions on the durations.
  @assert typeof(ruleDuration) <: Year
  @assert typeof(guarantee) <: Year

  # Ensure the output folder is ready.
  if length(saveFolder) > 0 && ! isdir(saveFolder)
    mkdir(saveFolder)
  end

  ### Generate parameters.
  hasCurrentRuleCurve = ! isempty(currentRuleCurve)

  generateDefaultLengthScenarios = in(:minimumRuleCurve, solvers) # Default length scenarios (i.e. guarantee).
  generateExtraLengthScenarios = in(:safeMinimumRuleCurve, solvers) # Longer scenarios (i.e. ruleDuration + guarantee).
  scenariosDefaultLength = Int(guarantee)
  scenariosExtraLength = Int(ruleDuration) + Int(guarantee)

  ruleLengthPeriods = Int(floor(Dates.days(ruleDuration) / Dates.days(getPeriod(reservoir))))
  guaranteePeriods = Int(floor(Dates.days(guarantee) / Dates.days(getPeriod(reservoir))))

  generateCIScenarios = length(ciAlgorithmsLevels) > 0
  ciAlgorithms = [algo[1] for algo in ciAlgorithmsLevels]

  # TODO: COMPLETE REFACTORING OF THE WHOLE FUNCTION! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  # A complete data structure of available solvers, presented with a uniform syntax:
  #     (reservoir, verbosity) -> f(...)::ReservoirSolution
  # TODO: Not yet useable throughout this script (just for cross-validation): for the solving process, needs information for
  # TODO:   the solution name (could just currify the functions for cross-validation and let that string as a parameter)
  # TODO:       (reservoir, verbosity, name) -> f(...)::ReservoirSolution
  # TODO:   however, this is still not enough, as minimumRuleCurve and safeMinimumRuleCurve do not take their scenarios from
  # TODO:   the same place (scgen_default, scgen_extra). Would need a more complex data structure to hold the solvers?
  # solversF = Dict{Symbol, Function}
  # for solver in solvers
  #   if solver == :minimumRuleCurve
  #     solversF[solver] = (sc, verbosity) -> minimumRuleCurve(sc, "", outputLength=ruleLengthPeriods, verbosity=verbosity)
  #   elseif solver == :safeMinimumRuleCurve
  #     solversF[solver] = (sc, verbosity) -> safeMinimumRuleCurve(sc, ruleLengthPeriods, guaranteePeriods, "", verbosity=verbosity)
  #   else
  #     error("Unrecognised solver: " * string(solver))
  #   end
  # end
  # TODO: COMPLETE REFACTORING OF THE WHOLE FUNCTION! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

  ### Handle scenario generation.
  if verbosity >= 0; println("[=*=] Scenario generation: started... "); end
  cis = Dict{Tuple{Symbol, Float64}, Tuple{Reservoir, Reservoir, Reservoir}}()
  scgen_default = Dict{Symbol, Reservoir}()
  scgen_extra = Dict{Symbol, Reservoir}()

  # Confidence interval computations if needed.
  for algo in ciAlgorithms
    for level in ciAlgorithmsLevels[algo]
      cis[(algo, level)] = computeCI(reservoir, level, ci_period=ciPeriod, method=algo)
    end
  end

  # Actual scenario generation.
  if generateDefaultLengthScenarios
    for algo in scenarioGenerationAlgorithms
      scgen_default[algo] = scenarioGeneration(reservoir, algo, scenariosDefaultLength)
    end
  end

  if generateExtraLengthScenarios
    for algo in scenarioGenerationAlgorithms
      scgen_extra[algo] = scenarioGeneration(reservoir, algo, scenariosExtraLength)
    end
  end

  if generateCIScenarios
    for algo in ciAlgorithms
      for level in ciAlgorithmsLevels[algo]
        s = symbol("ci_", algo, "_", level)
        if generateDefaultLengthScenarios
          scgen_default[s] = scenarioGeneration(cis[(algo, level)][2], :Duplicate, scenariosDefaultLength)
        end
        if generateExtraLengthScenarios
          scgen_extra[s] = scenarioGeneration(cis[(algo, level)][2], :Duplicate, scenariosExtraLength)
        end
      end
    end
  end

  if verbosity >= 0; println("[=*=] Scenario generation: done! "); end

  ### Handle solving.
  if verbosity >= 0; println("[=*=] Solving: started... "); end
  solutions = Dict{Tuple{Symbol, Symbol}, ReservoirSolution}()
  if hasCurrentRuleCurve
    solutions[(:currentRuleCurve, :currentRuleCurve)] = wrapExistingSolution(reservoir, currentRuleCurve)
  end

  for solver in solvers
    if verbosity >= 0; println("[=*=]   Solving with the algorithm: " * string(solver)); end

    if solver == :minimumRuleCurve
      for (lbl, sc) in scgen_default
        if ! in((lbl, solver), disabledCombinations)
          if verbosity >= 0; println("[=*=]     Solving with the data set: " * string(lbl)); end
          sol_lbl = "stochastic, default solver, " * string(guarantee) * ", " * string(lbl)
          solutions[(solver, lbl)] = minimumRuleCurve(sc, sol_lbl, outputLength=ruleLengthPeriods, useOutputs=getOutputIds(sc), verbosity=verbosity)
        end
      end
    elseif solver == :safeMinimumRuleCurve
      for (lbl, sc) in scgen_extra
        if ! in((lbl, solver), disabledCombinations)
          if verbosity >= 0; println("[=*=]     Solving with the data set: " * string(lbl)); end
          sol_lbl = "stochastic, safety solver, " * string(guarantee) * " + " * string(ruleDuration) * ", " * string(lbl)
          solutions[(solver, lbl)] = safeMinimumRuleCurve(sc, ruleLengthPeriods, guaranteePeriods, sol_lbl, verbosity=verbosity)
        end
      end
    else
      error("Unrecognised solver: " * string(solver))
    end
  end

  if verbosity >= 0; println("[=*=] Solving: done! "); end

  ## Check the consistency of the solutions.
  # All must have the same number of timesteps.
  sol_n_ts = -1 # Will be initialised to a positive value when the first feasible solution is found.
  sol_ts_pb = Dict{Any, Int}()
  for (desc, solution) in solutions
    if solution.feasible
      l = countTimeSteps(solution)
      if sol_n_ts < 0
        sol_n_ts = l
      elseif sol_n_ts != l
        sol_ts_pb[desc] = l
      end
    end
  end
  if sol_n_ts > 0 && length(sol_ts_pb) > 0
    txt_sol = "A solution has"
    if length(sol_ts_pb) > 1
      txt_sol = string(length(sol_ts_pb)) * " solutions have"
    end
    error(txt_sol * " an inconsistent number of time steps: " * string(sol_ts_pb))
  end

  # Other tests to perform?

  if verbosity >= 1; println("[=*=] Solving: all solutions are consistent! "); end

  ### Handle evaluation.
  if verbosity >= 0; println("[=*=] Evaluating: started... "); end

  evaluations = Dict{Tuple{Tuple{Symbol, Symbol}, Symbol}, ReservoirSolution}() # ((solver, scenarioGenerator), evaluator) -> evaluation
  # Associates (solutionAlgorithm, evaluator) to the result of an evaluation.
  # A solutionAlgorithm is a pair (solver, scenarioGenerator).
  for evaluator in evaluators
    if verbosity >= 0; println("[=*=]   Evaluating with the algorithm: " * string(evaluator)); end

    for (desc, solution) in solutions
      if ! solution.feasible
        continue
      end

      solver = desc[1]
      sc_gen = desc[2]

      if verbosity >= 0; println("[=*=]     Evaluating the solution: solver " * string(solver) * ", scenario generation " * string(sc_gen)); end

      if evaluator == :feasibilityProbabilitySyntheticScenarios
        evaluations[(desc, evaluator)] = feasibilityProbability(solution, numberScenarios=1000, scenarioDuration=guarantee, useOutputs=getOutputIds(getReservoir(solution)), verbosity=verbosity)
      elseif evaluator == :feasibilityProbabilityAllDatasets
        if verbosity >= 0; println("[=*=]     Evaluating with default length data (" * string(guarantee) * ")."); end
        for (lbl, sc) in scgen_default
         evaluations[(desc, symbol(evaluator, "_default_", lbl))] = feasibilityProbability(solution, scenarios=sc, useOutputs=getOutputIds(getReservoir(solution)), verbosity=verbosity)
        end
        if verbosity >= 0; println("[=*=]     Evaluating with extra length data (" * string(ruleDuration + guarantee) * ")."); end
        for (lbl, sc) in scgen_extra
          evaluations[(desc, symbol(evaluator, "_extra_", lbl))] = feasibilityProbability(solution, scenarios=sc, useOutputs=getOutputIds(getReservoir(solution)), verbosity=verbosity)
        end
      elseif evaluator == :feasibilityProbabilityClimateChange
        if verbosity >= 0; println("[=*=]     Evaluating with climate change scenarios."); end
        if isempty(climateChangeScenarios); error("No climate change scenarios given, while the script was asked to test against them!"); end;
        for (lbl, sc) in climateChangeScenarios
          evaluations[(desc, symbol(evaluator, "_", lbl))] = feasibilityProbability(solution, scenarios=sc, useOutputs=getOutputIds(getReservoir(solution)), verbosity=verbosity)
        end
      elseif evaluator == :crossValidatedFeasibilityProbability
        error("Cross-validation is not a solution evaluator! ")
      elseif evaluator == :worstCaseFunctions
        evaluations[(desc, evaluator)] = worstCaseFunctions(solution, useOutputs=getOutputIds(getReservoir(solution)), verbosity=verbosity)
      elseif evaluator == :maximumDeviation
        inflow = vec(getSolution(worstCaseFunctions(solutions[(:currentRuleCurve, :currentRuleCurve)], useOutputs=getOutputIds(getReservoir(solution)), verbosity=verbosity)))
        evaluations[(desc, evaluator)] = maximumDeviation(solution, inflow, verbosity=verbosity)
      elseif evaluator == :hydropower
        if solver == :minimumRuleCurve
          price_scenario = [sin(2 * π / 20 * t) for t in 1:scenariosDefaultLength]
        elseif solver == :safeMinimumRuleCurve
          price_scenario = [sin(2 * π / 30 * t) for t in 1:scenariosExtraLength]
        else
          error("Unrecognised solver for generating price scenarios: " * string(solver))
        end
        evaluations[(desc, evaluator)] = hydropower(solution, useOutputs=getOutputIds(getReservoir(solution)), electricityPrices=price_scenario)
      else
        error("Unrecognised evaluator: " * string(evaluator))
      end
    end
  end

  globalEvaluations = Dict{Tuple{Symbol, Symbol, Symbol}, ReservoirSolution}() # (solver, sc_gen, globalEvaluator) -> evaluation result
  for globalEvaluator in globalEvaluators
    if verbosity >= 0; println("[=*=]   Globally evaluating with the algorithm: " * string(globalEvaluator)); end

    for solver in solvers
      if verbosity >= 0; println("[=*=]     Evaluating solver with " * string(globalEvaluator) * ": " * string(solver)); end

      if globalEvaluator == :crossValidatedFeasibilityProbability
        if solver == :minimumRuleCurve
          for sc_gen in scenarioGenerationAlgorithms
            if in((sc_gen, solver), disabledCombinations)
              continue
            end

            if verbosity >= 0; println("[=*=]     Evaluating solver with " * string(globalEvaluator) * ": " * string(solver) * ", scenario generation " * string(sc_gen)); end

            f = function (sc, verbosity)
              scgen = scenarioGeneration(sc, sc_gen, scenariosDefaultLength)
              return minimumRuleCurve(scgen, "", useOutputs=getOutputIds(sc), outputLength=ruleLengthPeriods, verbosity=verbosity)
            end
            globalEvaluations[(solver, sc_gen, globalEvaluator)] = crossValidatedFeasibilityProbability(reservoir, f, verbosity=verbosity, folds=10)
          end
          for algo in ciAlgorithms
            for level in ciAlgorithmsLevels[algo]
              sc_gen = symbol("ci_", algo, "_", level)

              if in((sc_gen, solver), disabledCombinations)
                continue
              end

              if verbosity >= 0; println("[=*=]     Evaluating solver with " * string(globalEvaluator) * ": " * string(solver) * ", scenario generation " * string(sc_gen)); end

              f = function (sc, verbosity)
                cis = computeCI(sc, level, ci_period=ciPeriod, method=algo)[2]
                scgen = scenarioGeneration(cis, :Duplicate, scenariosDefaultLength)
                return minimumRuleCurve(scgen, "", useOutputs=getOutputIds(sc), outputLength=ruleLengthPeriods, verbosity=verbosity)
              end
              globalEvaluations[(solver, sc_gen, globalEvaluator)] = crossValidatedFeasibilityProbability(reservoir, f, verbosity=verbosity, folds=10)
            end
          end
        elseif solver == :safeMinimumRuleCurve
          for sc_gen in scenarioGenerationAlgorithms
            if in((sc_gen, solver), disabledCombinations)
              continue
            end

            if verbosity >= 0; println("[=*=]     Evaluating solver with " * string(globalEvaluator) * ": " * string(solver) * ", scenario generation " * string(sc_gen)); end

            f = function (sc, verbosity)
              scgen = scenarioGeneration(sc, sc_gen, scenariosExtraLength)
              return safeMinimumRuleCurve(scgen, ruleLengthPeriods, guaranteePeriods, "", verbosity=verbosity)
            end
            globalEvaluations[(solver, sc_gen, globalEvaluator)] = crossValidatedFeasibilityProbability(reservoir, f, verbosity=verbosity, folds=10)
          end
          for algo in ciAlgorithms
            for level in ciAlgorithmsLevels[algo]
              sc_gen = symbol("ci_", algo, "_", level)

              if in((sc_gen, solver), disabledCombinations)
                continue
              end

              if verbosity >= 0; println("[=*=]     Evaluating solver with " * string(globalEvaluator) * ": " * string(solver) * ", scenario generation " * string(sc_gen)); end

              f = function (sc, verbosity)
                cis = computeCI(sc, level, ci_period=ciPeriod, method=algo)[2]
                scgen = scenarioGeneration(cis, :Duplicate, scenariosExtraLength)
                return safeMinimumRuleCurve(scgen, ruleLengthPeriods, guaranteePeriods, "", verbosity=verbosity)
              end
              globalEvaluations[(solver, sc_gen, globalEvaluator)] = crossValidatedFeasibilityProbability(reservoir, f, verbosity=verbosity, folds=10)
            end
          end
        else
          error("Unrecognised solver: " * string(solver))
        end
      else
        error("Unrecognised global evaluator: " * string(solver))
      end
    end
  end

  if verbosity >= 0; println("[=*=] Evaluating: done! "); end

  ### Finally, display report, generate figures, and save what should be saved.
  if verbosity >= 1
    for desc in keys(solutions)
      println("**************************************************************************************************************")
      println("** Solver: " * string(desc[1]))
      println("** Scenario generation: " * string(desc[2]))
      println("**")
      for evaluator in evaluators
        # The generated keys for those solvers are not detected by the following condition!
        if evaluator == :feasibilityProbabilityAllDatasets
          for (eval_desc, evaluation) in evaluations
            if desc == eval_desc[1] && startswith(string(eval_desc[2]), "feasibilityProbabilityAllDatasets")
              ev = evaluation.solution
              first_part_length = length("feasibilityProbabilityAllDatasets_")
              dataset = string(eval_desc[2])[first_part_length + 1:end]
              println("** Feasibility analysis on data set " * dataset * ": probability of being infeasible: " * string(ev.probabilityInfeasible * 100) * "%")
              println("** Feasibility analysis on data set " * dataset * ": probability of being infeasible due to water surplus: " * string(ev.probabilityInfeasibleSurplus * 100) * "%")
              println("** Feasibility analysis on data set " * dataset * ": probability of being infeasible due to water shortage: " * string(ev.probabilityInfeasibleShortage * 100) * "%")
              println("** Feasibility analysis on data set " * dataset * ": average water shortage: " * string(ev.averageShortage) * " 10^6 m^3")
              println("** Feasibility analysis on data set " * dataset * ": average shortage index: " * string(ev.averageShortageIndex))
              println("** Feasibility analysis on data set " * dataset * ": average shortage vulnerability index: " * string(ev.averageShortageVulnerability))
            end
          end
        end
        if evaluator == :feasibilityProbabilityClimateChange
          for (eval_desc, evaluation) in evaluations
            if desc == eval_desc[1] && startswith(string(eval_desc[2]), "feasibilityProbabilityClimateChange")
              ev = evaluation.solution
              first_part_length = length("feasibilityProbabilityClimateChange_")
              dataset = string(eval_desc[2])[first_part_length + 1:end]
              println("** Feasibility analysis on the climate change scenario " * dataset * ": probability of being infeasible: " * string(ev.probabilityInfeasible * 100) * "%")
              println("** Feasibility analysis on the climate change scenario " * dataset * ": probability of being infeasible due to water surplus: " * string(ev.probabilityInfeasibleSurplus * 100) * "%")
              println("** Feasibility analysis on the climate change scenario " * dataset * ": probability of being infeasible due to water shortage: " * string(ev.probabilityInfeasibleShortage * 100) * "%")
              println("** Feasibility analysis on the climate change scenario " * dataset * ": average water shortage: " * string(ev.averageShortage) * " 10^6 m^3")
              println("** Feasibility analysis on the climate change scenario " * dataset * ": average shortage index: " * string(ev.averageShortageIndex))
              println("** Feasibility analysis on the climate change scenario " * dataset * ": average shortage vulnerability index: " * string(ev.averageShortageVulnerability))
            end
          end
        end

        if haskey(evaluations, (desc, evaluator)) # Not all evaluators could work previously (e.g. cross-validation needs more scenarios than available with CIs, with more parameters than a sheer symbol).
          ev = evaluations[(desc, evaluator)].solution
          if evaluator == :feasibilityProbabilitySyntheticScenarios
            println("** Monte-Carlo feasibility analysis: probability of being infeasible: " * string(ev.probabilityInfeasible * 100) * "%")
            println("** Monte-Carlo feasibility analysis: probability of being infeasible due to water surplus: " * string(ev.probabilityInfeasibleSurplus * 100) * "%")
            println("** Monte-Carlo feasibility analysis: probability of being infeasible due to water shortage: " * string(ev.probabilityInfeasibleShortage * 100) * "%")
            println("** Monte-Carlo feasibility analysis: average water shortage: " * string(ev.averageShortage) * " 10^6 m^3")
            println("** Monte-Carlo feasibility analysis: average shortage index: " * string(ev.averageShortageIndex))
            println("** Monte-Carlo feasibility analysis: average shortage vulnerability index: " * string(ev.averageShortageVulnerability))
          elseif evaluator == :worstCaseFunctions
            println("** Inverse model analysis: annual volume below which infeasibility is guaranteed: " * string(ev.total) * " 10^6 m^3")
          elseif evaluator == :maximumDeviation
            println("** Maximum deviation: from the inverse solution for the current rule curve: " * string(ev.deviation) * " 10^6 m^3")
          else
            error("Unrecognised evaluator: " * string(evaluator))
          end
        end
      end
      println("**************************************************************************************************************")
    end

    for globalEvaluator in globalEvaluators
      if globalEvaluator == :crossValidatedFeasibilityProbability
        for (globalEvaluationDescription, globalEvaluation) in globalEvaluations
          solver = globalEvaluationDescription[1]
          sc_gen = globalEvaluationDescription[2]
          evaluator = globalEvaluationDescription[3]
          if evaluator == :crossValidatedFeasibilityProbability
            ev = globalEvaluation.solution

            println("**************************************************************************************************************")
            println("** Solver: ", solver)
            println("** Scenario generation: ", sc_gen)
            println("** ")
            println("** Cross-validation feasibility analysis: probability of being infeasible: " * string(ev.probabilityInfeasible * 100) * "%")
            println("** Cross-validation feasibility analysis: probability of being infeasible due to water surplus: " * string(ev.probabilityInfeasibleSurplus * 100) * "%")
            println("** Cross-validation feasibility analysis: probability of being infeasible due to water shortage: " * string(ev.probabilityInfeasibleShortage * 100) * "%")
            println("** Cross-validation feasibility analysis: probability of having no rule curve: " * string(ev.probabilityNoRuleCurve * 100) * "%")
            println("** Cross-validation feasibility analysis: average water shortage: " * string(ev.averageShortage) * " 10^6 m^3")
            println("** Cross-validation feasibility analysis: average shortage index: " * string(ev.averageShortageIndex))
            println("** Cross-validation feasibility analysis: average shortage vulnerability index: " * string(ev.averageShortageVulnerability))
            println("**************************************************************************************************************")
          end
        end
      end
    end
  end

  if saveIntermediate
    if verbosity >= 0; println("[=*=] Serialising: started... "); end

    folder_intermediate = saveFolder * "/intermediate/"
    if ! isdir(folder_intermediate)
      mkdir(folder_intermediate)
    end
    serialise_own(filename, o) = open(folder_intermediate * filename, "w") do stream; serialize(stream, o); end

    serialise_own("cis.jls", cis)
    serialise_own("scgen_default.jls", scgen_default)
    serialise_own("scgen_extra.jls", scgen_extra)
    serialise_own("solutions.jls", solutions)
    serialise_own("evaluations.jls", evaluations)
    serialise_own("global_evaluations.jls", globalEvaluations)

    if verbosity >= 0; println("[=*=] Serialising: done! "); end
  end

  if savePlots
    if verbosity >= 0; println("[=*=] Plotting: started... "); end

    folder_plots = saveFolder * "/plots/"
    if ! isdir(folder_plots)
      mkdir(folder_plots)
    end

    function exportPlot(plot::AbstractPlot, name::AbstractString)
      savefig(plot, folder_plots * name * ".png")
      savefig(plot, folder_plots * name * ".svg")
      nothing
    end

    function exportData(data::Array{Any, 2}, name::AbstractString)
      writedlm(folder_plots * name * ".tsv", data, '\t')
      nothing
    end

    period = getPeriod(reservoir)

    ## Data treatment.
    # Confidence intervals. First plots for each algorithm, then a global plot.
    # Indices in the code: once the average [1], then for each pair upper/lower bounds [2] and [3].
    if ! isempty(currentRuleCurve) && length(ciAlgorithms) > 0
      # Plot each algorithm independently (if at least two algorithms), all confidence levels at once.
      if length(ciAlgorithms) > 1
        for algo in ciAlgorithms
          # All confidence levels at once.
          data = vcat(cis[(algo, ciAlgorithmsLevels[algo][1])][1], [let lcis = cis[(algo, level)]; [lcis[2], lcis[3]]; end for level in ciAlgorithmsLevels[algo]]...)
          plt, dat = plotRiverAlongReservoirs(data, period, 1:ruleLengthPeriods, getName(reservoir), exportPlotData=true)
          exportPlot(plt, "confidence_intervals_" * string(algo))
          exportData(dat, "confidence_intervals_" * string(algo))
        end
      end

      # Each confidence level in its own plot (if more than one level requested), for each algorithm.
      for algo in ciAlgorithms
        if length(ciAlgorithmsLevels[algo]) > 1
          for level in ciAlgorithmsLevels[algo]
            data = vcat(cis[(algo, ciAlgorithmsLevels[algo][1])][1], cis[(algo, level)][2], cis[(algo, level)][3])
            plt, dat = plotRiverAlongReservoirs(data, period, 1:ruleLengthPeriods, getName(reservoir), exportPlotData=true)
            exportPlot(plt, "confidence_intervals_" * string(algo) * "_cl" * string(level))
            exportData(dat, "confidence_intervals_" * string(algo) * "_cl" * string(level))
          end
        end
      end

      # TODO: plot for a given confidence level (irrespectively of algorithms).

      # Global plot for all algorithms and confidence levels at once.
      data = vcat([[cis[(algo, ciAlgorithmsLevels[algo][1])][1]; vcat([let lcis = cis[(algo, level)]; [lcis[2]; lcis[3]]; end for level in ciAlgorithmsLevels[algo]]...)] for algo in ciAlgorithms]...)
      plt, dat = plotRiverAlongReservoirs(data, period, 1:ruleLengthPeriods, getName(reservoir), exportPlotData=true)
      exportPlot(plt, "confidence_intervals")
      exportData(dat, "confidence_intervals")
    end

    # Skip scenario generation output.

    ## Solutions.
    # Existing target, alone.
    if ! isempty(currentRuleCurve)
      plt, dat = plotRule(currentRuleCurve, period, title="Current rule curve", exportPlotData=true)
      exportPlot(plt, "current_rule")
      exportData(dat, "current_rule")
    end

    # Complete set of scenarios when merging them and compare them to the solution.
    if in(:Merge, scenarioGenerationAlgorithms) && in(:minimumRuleCurve, solvers) && ! isempty(currentRuleCurve)
      merge_sol = solutions[(:minimumRuleCurve, :Merge)]
      colourings = [:YearlyAverage, :WetSeasonAverage, :DrySeasonAverage, :DriestMonth, :DriestThreeMonths, :DriestSixMonths]

      for colouring in colourings
        plt, dat = plotSolutionVsScenarios(merge_sol, currentRuleCurve, colouring=colouring, exportPlotData=true)
        exportPlot(plt, "solution_safety_stochastic_merge_splitperscenario_" * string(colouring))
        exportData(dat, "solution_safety_stochastic_merge_splitperscenario_" * string(colouring))
      end
    end

    # Solutions for the various solvers: one plot for each solver first, then one plot with each and every solution.
    if ! isempty(currentRuleCurve)
      for solver in solvers
        solver_solution = filter((k, v) -> k[1] == solver, solutions)
        plt, dat = plotSolutions(collect(values(solver_solution)), currentRuleCurve, exportPlotData=true)
        exportPlot(plt, "solutions_" * string(solver))
        exportData(dat, "solutions_" * string(solver))
      end

      plt, dat = plotSolutions(collect(values(solutions)), currentRuleCurve, exportPlotData=true)
      exportPlot(plt, "solutions")
      exportData(dat, "solutions")
    end

    ## Evaluation.
    # Nothing for solution evaluators.
    # Nothing for global evaluators.

    if verbosity >= 0; println("[=*=] Plotting: done! "); end
  end

  return cis, scgen_default, scgen_extra, solutions, evaluations
end
