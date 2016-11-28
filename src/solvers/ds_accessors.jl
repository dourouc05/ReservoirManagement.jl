["""
Getters for the solver data structures.
"""]



getReservoir(sol::ReservoirSolution) = sol.reservoir
isFeasible(sol::ReservoirSolution) = sol.feasible
getSolutionObject(sol::ReservoirSolution) = sol.solution
hasName(sol::ReservoirSolution) = isnull(sol.name)
getName(sol::ReservoirSolution) = get(sol.name)
getSolverIdentifier(sol::ReservoirSolution) = sol.solver
getSolverDescription(sol::ReservoirSolution) = sol.solverDescription
getSolverOptions(sol::ReservoirSolution) = sol.solverOptions
getSolverStatistics(sol::ReservoirSolution) = sol.solverStatistics

isInfeasible(sol::ReservoirSolution) = ! isFeasible(sol)
getNameOrElse(sol::ReservoirSolution) = get(sol.name, getVariant(sol))
getNameOrElse(sol::ReservoirSolution, default::AbstractString) = get(sol.name, default)

getPeriod(sol::ReservoirSolution) = getPeriod(getReservoir(sol))
getVariant(sol::ReservoirSolution) = getVariant(getReservoir(sol))
getSolution(sol::ReservoirSolution) = getSolution(getSolutionObject(sol)) # Dive directly into the data structures, as opposed to getSolutionObject.
countTimeSteps(sol::ReservoirSolution) = countTimeSteps(getSolutionObject(sol))


getSolution(sol::Solution) = sol.solution
function countTimeSteps(sol::Solution)
  error("countTimeSteps not implemented for " * string(typeof(sol)))
end
countTimeSteps(sol::EvaluationSolution) = 0 # This usually has no sense for evaluation, as the result is mainly one number. Can still be overriden. 
