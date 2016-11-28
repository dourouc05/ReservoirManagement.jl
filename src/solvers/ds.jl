["""
Defines generic solver facilities. Hence all data structures here must be abstract! Here, a "solver" is defined
as some piece of code that takes some reservoir in input and outputs some results based on this reservoir.

A solution is made up of four parts:

  * the main solution object, whose type is `ReservoirSolution`, and includes the three other ones.
    It remembers the input data and general information, plus solver-specific information.
  * a solution object, sub-type of `Solution`, storing the output from this solver.
  * an option object, sub-type of `Options`, remembering the options the solver was called with.
  * a statistics object, sub-type of `SolverStatistics`, which collects statistics about the solving process.

Those three solver-dependent objects follow the same hierarchy as the solvers, which goes as follows
(the first level take into account the differences in output signification and input characteristics):

  * The *rule curve solvers*, whose output is a rule curve, based on a series of scenarios.
      * *Minimum* rule curve:
          * `minimumRuleCurve`: computes the minimum level in the reservoir so it can fulfil all its functions
            for the optimisation time horizon (i.e. guarantees correct operations for decreasing durations)
          * `safeMinimumRuleCurve`: computes the minimum level in the reservoir so it can fulfil all its functions
            for any given time horizon (i.e. guarantees correct operations for any considered period with
            the same level of safety)
      * *Maximum* rule curve:
          * `safeMaximumRuleCurve`: computes the maximum level in the reservoir so it can withstand spate events
            for any given time horizon (i.e. guarantees correct operations for any considered period with
            the same level of safety)
  * The *event solvers*, whose output is the behaviour for one specific event.
      * `spateMaximumLevel`: computes the maximum level in the reservoir so it can withstand the given spate event.
        This is not a rule curve, as it is valid only for the specific event.
        TODO!
  * The *inverse solvers*, whose output is a scenario, evaluating a rule curve or any other decision.
      * `worstCaseFunctions`: computes the worst-case inflow scenario so the reservoir is hardly able to fulfil
        its functions.
  * The *evaluation solvers*, whose output is a quantity evaluating a rule curve.
      * `feasibilityProbability`: computes a probability of failure for a given solution by generating new scenarios
        with a statistical model.
      * `purposeShortage`: determines whether the infeasibility of a solution is due to a shortage of water
        (i.e. impossible to meet all purposes) or a surplus thereof (i.e. a flood issue).
      * `crossValidatedFeasibilityProbability`: computes a probability of failure for an algorithm, i.e. reoptimises
        a solution for various subsets of the data. (Useful when data is not profuse.)
      * `maximumDeviation`, whose output is a maximum allowable deviation with respect to a given inflow scenario.
      * `hydropower`, whose output is a hydropower potential for a given solution.

New solver
==========

To write a new solver:

  * Prepare a subtype of `Options` to store the various options
  * Prepare a subtype of `Solution` to store all interesting properties of the solution (including by-products)
  * Prepare a subtype of `Statistics` to store all interesting properties of the solving process (such as required time, number of iterations)
  * Write a solver function, whose input starts with a reservoir or a solutio, and may be completed with other parameters.
    Its output must be a `ReservoirSolution` using the previously defined data structures.

The data structures should be included within the base nodes in this file:

  * Taking a reservoir as input:
    * `RuleCurve*` for those optimising a situation to get a rule curve (either as a mandatory level,
      or bounds on the level)
      * `MinimumRuleCurve*` more specifically for those imposing a minimum level for some purpose
      * `MaximumRuleCurve*` more specifically for those imposing a maximum level for some purpose
    * `Event*` for those optimising a situation to get a rule curve that best responds to some event (such as flooding)
  * Taking a solution as input:
    * `Inverse*` for the inverse solvers, working from a solution to the minimum inflow to reach this minimum rule curve, exactly
    * `Evaluation*` for the evaluations of a solution, be it in terms of probability or other assessments
"""]



## Generic data structures.

"""
Generic abstraction of solver options. All options given to a solver with a given drainage basin
to specify completely its behaviour.
"""
abstract Options

"""
Generic abstraction of solver solution. All subtypes contain at least a `solution` field that contains the solution,
but each solver may consider other variables for inclusion.
"""
abstract Solution

"""
Time for a solver to get a solution, expressed in milliseconds. All have a `total` field expressing the total time the solver took.
"""
abstract SolverStatistics

"""
The *context* of a solution for a reservoir after being solved by some solver.

Available fields:

  * `reservoir`: the reservoir that was used by the solver
  * `solution`: the actual solution given by the solver (the exact data structure depends on the solver; it must be a subtype of `Solution`)
  * `solver`: a symbol indicating what solver was used; this field can be used to perform solver identification
  * `solverDescription`: a brief textual description of the solver
  * `solverOptions`: further specifications on the solver that is used (the exact data structure depends on the solver; it must be a subtype of `Options`)
  * `solverStatistics`: a description of the runtime aspects of the solving process, such as timings (the exact data structure depends on the solver; it must be a subtype of `SolverStatistics`)
"""
immutable ReservoirSolution
  reservoir::Reservoir
  feasible::Bool # TODO: does it make sense here for all solvers? NO!
  solution::Solution
  name::Nullable{ASCIIString}

  solver::Symbol
  solverDescription::AbstractString
  solverOptions::Options
  solverStatistics::SolverStatistics
end

Base.isempty(r::ReservoirSolution) = isempty(r.solution)



## Rule curve.

"""
Generic abstraction for the options of a solver optimising a rule curve.

  * objective::Symbol
"""
abstract RuleCurveOptions <: Options

"""
A more specific kind of solution that takes the form of a rule curve over time.
"""
abstract RuleCurveSolution <: Solution

immutable EmptyRuleCurveSolution <: RuleCurveSolution end
Base.isempty(s::EmptyRuleCurveSolution) = true
Base.isempty(s::RuleCurveSolution) = false

"""
A more specific kind of solver statistics for rule curves.
"""
abstract RuleCurveSolverStatistics <: SolverStatistics



# Minimum rule curve.

"""
Generic abstraction for the options of a solver optimising a minimum rule curve.

  * scenarioProbabilities::Array{Float64, 1}
  * objective::Symbol
  * imposeDiverted::Bool
"""
abstract MinimumRuleCurveOptions <: RuleCurveOptions

"""
A more specific kind of solution that takes the form of a minimum rule curve over time.
"""
abstract MinimumRuleCurveSolution <: RuleCurveSolution

immutable EmptyMinimumRuleCurveSolution <: MinimumRuleCurveSolution end
Base.isempty(s::EmptyMinimumRuleCurveSolution) = true
Base.isempty(s::MinimumRuleCurveSolution) = false

"""
A more specific kind of solver statistics for minimum rule curves.
"""
abstract MinimumRuleCurveSolverStatistics <: RuleCurveSolverStatistics



# Maximum rule curve.

"""
Generic abstraction for the options of a solver optimising a maximum rule curve.

  * objective::Symbol
"""
abstract MaximumRuleCurveOptions <: RuleCurveOptions

"""
A more specific kind of solution that takes the form of a maximum rule curve over time.
"""
abstract MaximumRuleCurveSolution <: RuleCurveSolution

immutable EmptyMaximumRuleCurveSolution <: MaximumRuleCurveSolution end
Base.isempty(s::EmptyMaximumRuleCurveSolution) = true
Base.isempty(s::MaximumRuleCurveSolution) = false

"""
A more specific kind of solver statistics for maximum rule curves.
"""
abstract MaximumRuleCurveSolverStatistics <: RuleCurveSolverStatistics



## Event.

"""
Generic abstraction for the options of a solver optimising a response to a given event.

  * objective::Symbol
"""
abstract EventOptions <: Options

"""
A more specific kind of solution that takes the form of a response to a given event.
"""
abstract EventSolution <: Solution

immutable EmptyEventSolution <: EventSolution end
Base.isempty(s::EmptyEventSolution) = true
Base.isempty(s::EventSolution) = false

"""
A more specific kind of solver statistics for event responses.
"""
abstract EventSolverStatistics <: SolverStatistics



## Inverse.

"""
Generic abstraction for the options of an inverse solver.
"""
abstract InverseOptions <: Options

"""
A more specific kind of solution that takes the form of an input to other models (i.e. after an inverse model).
"""
abstract InverseSolution <: Solution

immutable EmptyInverseSolution <: InverseSolution end
Base.isempty(s::EmptyInverseSolution) = true
Base.isempty(s::InverseSolution) = false

"""
A more specific kind of solver statistics for an inverse model.
"""
abstract InverseSolverStatistics <: SolverStatistics



## Evaluation.

"""
Generic abstraction for the options of an evaluation.
"""
abstract EvaluationOptions <: Options

"""
A more specific kind of solution that evaluates another one.
"""
abstract EvaluationSolution <: Solution

immutable EmptyEvaluationSolution <: EvaluationSolution end
Base.isempty(s::EmptyEvaluationSolution) = true
Base.isempty(s::EvaluationSolution) = false

"""
A more specific kind of solver statistics for an evaluation.
"""
abstract EvaluationSolverStatistics <: SolverStatistics
