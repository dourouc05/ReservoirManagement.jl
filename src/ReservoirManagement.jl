module ReservoirManagement

  using TimeSeries
  using Base.Dates
  using DataStructures
  using SIUnits
  using SIUnits.ShortUnits

  using JuMP
  using Gurobi

  using StatsBase
  using HypothesisTests
  using Distributions
  using ExtremeValueDistributions # Strange dependency (not published in METADATA), try to reimplement somewhere
  # (and make it fully compatible with Julia 0.4, without warnings).

  using DataFrames
  using Colors
  using Plots
  using StatPlots

  # The generalised distributions are defined in both Distributions and ExtremeValueDistributions, but those are identical.
  GeneralizedExtremeValue = Distributions.GeneralizedExtremeValue
  GeneralizedPareto = Distributions.GeneralizedPareto

  # Conflict between standard Dates and SIUnits.
  Second = Dates.Second

  # To allow overriding those functions.
  import Base: copy, isempty, start, next, done, eltype, length

  # Order of inclusion is important!
  include("utils.jl")
  include("reservoir/_ds.jl")
  include("read.jl")
  include("stats/_stats.jl")
  include("scenariosgeneration.jl")
  include("solvers/_solvers.jl")
  include("plots_dataframe.jl")
  include("plots_plots.jl")
  include("report.jl")

  # Export things that should be exported.
  export
    ## utils.jl
    Height, Volume, Discharge, Power,
    SeasonsSeparator, HydrologicalWesternEurope, seasons, isSeasonKnown, ensureSeasonKnown, begins, ends, length, seasonBoundariesGuarantee, separate,
    NonLeapYearIterator,
    ArrayAsIterableReservoirScenarios, eachscenario,
    ArrayAsIterableReservoirScenario, totalInflow,
    toLength,

    ## reservoir/
    # bathymetry.jl
    ReservoirBathymetry, NullReservoirBathymetry, computeVolume, volume, computeLevel, level,
    WhiteBoxReservoirBathymetry, iswhitebox, isblackbox,
    PolynomialReservoirBathymetry, ispolynomial, isPolynomial, coefficient,
    LinearReservoirBathymetry, getIntercept, getSlope, islinear, isLinear, intercept, slope,
    QuadraticReservoirBathymetry, isquadratic, isQuadratic,
    PiecewiseLinearReservoirBathymetry, ispiecewiselinear, isPiecewiseLinear, isPWL,
    SquareRootReservoirBathymetry, issquareroot, isSquareRoot,

    # hydropower.jl
    HydropowerUnit,
    getId, getPurposeId, getDamOutputId, getMaximumDischarge, getMaximumPower, getEfficiency,
    id, purposeId, damOutputId, maxDischarge, maxPower, efficiency,

    # purpose.jl
    Purpose, getId, isDeterministic, id,
    WaterWithdrawalPurpose, isWaterWithdrawal,
    DeterministicWaterWithdrawalPurpose, getNeed, need, DeterministicDrinkingWater, DeterministicEnvironmentalFlow,
    DeterministicHydropower, getHydropowerUnit, isHydropower, hydropower,
    DeterministicPurposes,

    # output.jl
    DamOutput, getId, getCondition, id, condition, hasCondition, isConditional, 
    ConveyanceDamOutput, hasConveyanceFactor, getConveyanceFactor, conveyance,
    NoCondition, isempty,
    ConstantDamOutput, getDischarge, isConstant, discharge,
    HydropowerDamOutput, getId, getDischarge, isHydropower, getHydropowerUnit, hydropower,
    ConveyanceSpillway, ConveyancePipe,
    MinimumReservoirLevelCondition, getMinimumLevel, minLevel,
    ConstantDamOutputs,

    # river.jl
    River, getName, getScenarios, getScenarioTS, getScenario,
    countScenarios, countTimeSteps, getPeriod, getTimePoints,
    name, scenarios, scenario, period,

    # river_natural.jl
    NaturalRiver,

    # river_diverted.jl
    DivertedRiver,
    getEnvironmentalFlow, getMaximumFlow, getMaximumAllowableFlow,
    environmentalFlow, maxFlow, maxAllowableFlow,

    # reservoir.jl
    Reservoir, getName, getVariant, getMinimumCapacity, getMaximumCapacity, getHydropowerUnits,
    hasHydropower, getHydropowerIds, getPenstockHydropower,
    getPurposes, getOutputs, getBathymetry, hasPurpose, getPurpose, getPurposeIds, getOutputIds, hasOutput, getOutput,
    name, variant, minCapacity, maxCapacity, penstockHydropower, hydropowerUnits, outputs, outputsId, output, bathymetry, purpose,
    countScenarios, countTimeSteps, getPeriod, getTimePoints,
    rivers, naturalRivers, divertedRivers, getNaturalRivers, getDivertedRivers, getRivers, getNaturalRiver, getDivertedRiver, hasRiver, getRiver,
    getDivertedRiverNames, getNaturalRiverNames, getRiverNames,
    totalInflow, removeScenarios, keepTimeSteps,
    NaturalRiverIterator, eachnaturalriver,
    DivertedRiverIterator, eachdivertedriver,
    RiverIterator, eachriver,
    ScenarioReservoirIterator, eachscenario,

    # basin.jl
    # Not for now: nothing uses it. For future multireservoir implementation.

    # bathymetry_fwd.jl
    level,

    # output_fwd.jl
    getMinimumVolume, minVolume,

    # misc.jl
    countScenarios, countTimeSteps,

    ## read.jl
    makeScenarios,

    ## stats/
    # base.jl
    computeCI,
    RiverScenarioIterator, eachscenario,
    generateScenarios,
    ReservoirScenarioIterator, eachscenario,

    # utils.jl
    splitIntoSets,

    # adam2015.jl
    printReport,

    # tstudent.jl
    # Nothing, no user-visible interface (only through the functions of base.jl).

    ## scenariosgeneration.jl
    scenarioDuplication, scenarioConcatenation, scenarioMerging, scenarioMixing,
    addScenarioGenerationAlgorithm!,
    scenarioGeneration, scenarioDuplication, scenarioConcatenation, scenarioConcatenation, scenarioMerging, scenarioMixing,
    scenarioShift,

    ## solvers/
    # ds.jl
    Options, Solution, SolverStatistics, ReservoirSolution,
    RuleCurveOptions, RuleCurveSolution, EmptyRuleCurveSolution, RuleCurveSolverStatistics,
    MinimumRuleCurveOptions, MinimumRuleCurveSolution, EmptyMinimumRuleCurveSolution, MinimumRuleCurveSolverStatistics,
    MaximumRuleCurveOptions, MaximumRuleCurveSolution, EmptyMaximumRuleCurveSolution, MaximumRuleCurveSolverStatistics,
    EventOptions, EventSolution, EmptyEventSolution, EventSolverStatistics,
    InverseOptions, InverseSolution, EmptyInverseSolution, InverseSolverStatistics,
    EvaluationOptions, EvaluationSolution, EmptyEvaluationSolution, EmptyEvaluationSolution, EvaluationSolverStatistics,

    # ds_accessors.jl
    getReservoir, isFeasible, getSolutionObject, hasName, getName, getSolverIdentifier, getSolverDescription, getSolverOptions, getSolverStatistics,
    isInfeasible, getNameOrElse, getNameOrElse,
    getPeriod, getVariant,
    getSolution, countTimeSteps,

    # none.jl
    NoSolverOptions, NoSolverSolution, NoSolverStatistics, wrapExistingSolution,

    # _helpers_inflow.jl
    # _helpers_output.jl
    # _helpers_bathymetry.jl
    # _helpers_purpose.jl
    # For now, helpers API not old enough to be exposed to the public.

    # evaluate_crossvalidate.jl
    FeasibilityCrossValidationOptions, FeasibilityCrossValidationSolution, FeasibilityCrossValidationSolverStatistics,
    crossValidatedFeasibilityProbability,

    # evaluate_deviation.jl
    DeviationOptions, DeviationSolution, DeviationSolverStatistics,
    maximumDeviation,

    # evaluate_mc_feasibility.jl
    FeasibilityMonteCarloOptions, FeasibilityMonteCarloSolution, FeasibilityMonteCarloSolverStatistics,
    feasibilityProbability,

    # evaluate_purposeshortage.jl
    ShortageOptions, ShortageSolution, ShortageSolverStatistics,
    computePurposeShortageIndex,

    # evaluate_hydropower.jl
    # Nothing for now, far from working state.

    # inverse_functions.jl
    FunctionsInverseOptions, FunctionsInverseSolution, FunctionsInverseSolverStatistics,
    worstCaseFunctions,

    # rulecurve_min_default.jl
    DefaultMinimumRuleCurveOptions, DefaultMinimumRuleCurveSolution, DefaultMinimumRuleCurveSolverStatistics,
    minimumRuleCurve,

    # rulecurve_min_safe.jl
    SafetyMinimumRuleCurveOptions, SafetyMinimumRuleCurveSolution, SafetyMinimumRuleCurveSolverStatistics,
    safeMinimumRuleCurve,

    ## plots_
    # dataframes.jl
    getDataFrameForRiver, getDataFrameForSolutions,

    # plots.jl
    plotRiverAlongReservoirs, plotSolutionVsScenarios, plotSolutions, plotRule,

    ## report.jl
    report
end
