["""
Gives a minimum rule curve with ensured safety guarantees.
"""]



immutable SafetyMinimumRuleCurveOptions <: MinimumRuleCurveOptions
  scenarioProbabilities::Array{Float64, 1}
  objective::Symbol
  imposeDiverted::Bool
  useOutputs::Array{Symbol, 1}

  outputLength::Int
  safetyGuarantee::Int
end

"""
The safety solver solution only has a rule curve as output.
"""
immutable SafetyMinimumRuleCurveSolution <: MinimumRuleCurveSolution
  solution::Array{Float64, 1}
end

countTimeSteps(sol::SafetyMinimumRuleCurveSolution) = length(sol.solution)

immutable SafetyMinimumRuleCurveSolverStatistics <: MinimumRuleCurveSolverStatistics
  total::Float64
  iterations::Array{DefaultMinimumRuleCurveSolverStatistics, 1}
end



"""
Ensure safety for a given period of time: the output is a rule that ensures feasibility of the solution for exactly
the asked period of time. (No points will have a longer safety guarantee, which would be the case with simple
optimisation procedure over a long period of time, like `minimumRuleCurve`.)
This is an application of model predictive control (MPC), more precisely receding horizon control
(http://www.springer.com/la/book/9781846280245).

The output is a rule over `outputLength`, with one point per time step as in `reservoir`. It must have at least
`outputLength` plus `safetyGuarantee` time steps of scenarios.

On the implementation side, this is equivalent to running the optimisation over a moving time horizon of the given
length, starting at various points within the first time step up to `outputLength`, with a horizon of always
`safetyGuarantee`.

Inputs:

  * `reservoir`: a situation to optimise
  * `outputLength`: a number of time steps for the output solution (with the same discretisation as `reservoir`)
  * `safetyGuarantee`: a number of time steps for the solution to remain feasible (with the same discretisation
                       as `reservoir`)
  * `filtering`: a function called at each iteration to filter the scenarios that will be used for this iteration,
                 mainly to remove scenarios that are too close to each other due to shifting the time horizon.
                 Prototype:
                     `(currentIter::Int, maxIter::Int, safetyGuarantee::Int, r::Reservoir, verbose::Bools) -> scenarios`
  * `verbosity`: level of output.
    * `-1`: no output
    * `0`: no output, except for solver's warnings
    * `1`: summary output
    * `2`: detailed output of this algorithm's details
    * `3`: detailed output of this algorithm's and its sub-algorithms' details
    * `4`: including solver output

TODO: replace outputLength and safetyGuarantee by periods of time.
TODO: check scenarios are long enough for the given lengths.
"""
function safeMinimumRuleCurve(reservoir::Reservoir, outputLength::Int, safetyGuarantee::Int, name::AbstractString="";
                              probabilities=zeros(Float64, 0), objective::Symbol=:Rule, imposeDiverted::Bool=true,
                              useOutputs::Array{Symbol, 1}=Symbol[],
                              filtering::Function=(currentIter::Int, maxIter::Int, safetyGuarantee::Int, r::Reservoir, verbose::Bool) -> r,
                              verbosity::Int=1)
  if verbosity >= 1
    println(" == ")
    println(" ==== Optimisation model for operational rules with MPC: ", getVariant(reservoir))
    println(" ==== MPC iterations to perform: ", outputLength)
  end

  sub_verbosity = (verbosity >= 3) ? verbosity : -1
  order_magnitude_iters = 10 ^ floor(Int, log10(outputLength))

  time_start = float(time_ns())

  output = zeros(Float64, outputLength)
  feasible = true
  statistics = Array(DefaultMinimumRuleCurveSolverStatistics, outputLength)
  for n = 1:outputLength
    if verbosity >= 1
      if verbosity >= 2 # Condition to choose the indentation.
        println(" ==== Iteration #", n, " out of ", outputLength)
      elseif n % order_magnitude_iters == 0
        println(" ====== Iteration #", n, " out of ", outputLength)
      end
    end

    # Compute the new solution.
    shifted_reservoir = filtering(n, outputLength, safetyGuarantee, scenarioShift(reservoir, firstKept=n, length=safetyGuarantee), sub_verbosity > 0)
    sol = minimumRuleCurve(shifted_reservoir,
                           probabilities=probabilities, objective=objective, imposeDiverted=imposeDiverted,
                           useOutputs=useOutputs,
                           verbosity=sub_verbosity) # TODO: prefix = "==" so that output is still legible (indented from this).
    target = sol.solution.solution
    statistics[n] = sol.solverStatistics

    # Check for feasibility and prepare corresponding output.
    if reduce(|, isnan(target)) # any(isnan(target))
      feasible = false
      output = zeros(size(output)) / 0.
      break
    end
    output[n] = target[1]

    if verbosity >= 3
      println(" ======== Iteration #", n, " provided point: ", output[n])
    end
  end

  time_done = float(time_ns())

  if verbosity >= 1
    println(" ==== Total storage: ", sum(output) * 1000000, " m^3 over ", outputLength * getPeriod(reservoir), " (safety guarantee: ", getPeriod(reservoir) * safetyGuarantee, ").")
    println(" ==== Total time for ", outputLength, " iterations: ", (time_done - time_start) / 1000000, "ms")
    println(" == ")
  end

  options = SafetyMinimumRuleCurveOptions(probabilities, objective, imposeDiverted, useOutputs, outputLength, safetyGuarantee)
  solution = SafetyMinimumRuleCurveSolution(output)
  time = SafetyMinimumRuleCurveSolverStatistics((time_done - time_start) / 1000000, statistics)
  solName = (name == "") ? Nullable{AbstractString}() : Nullable(name)
  return ReservoirSolution(reservoir, feasible, solution, solName, :SafeRuleCurveSolver, "safe rule curve from inflow", options, time)
end
