["""
Wraps a preexisting solution in the `ReservoirSolution` machinery.
"""]



immutable NoSolverOptions <: Options end
immutable NoSolverSolution <: Solution
  solution::Array{Float64, 1}
end
immutable NoSolverStatistics <: SolverStatistics end

Base.isempty(s::NoSolverSolution) = isempty(s.solution)
countTimeSteps(sol::NoSolverSolution) = length(sol.solution)

wrapExistingSolution(reservoir::Reservoir, solution::Array{Float64, 1}, name::AbstractString="Original rule curve") =
  ReservoirSolution(reservoir, true, NoSolverSolution(solution), Nullable(name), :NoSolver, "", NoSolverOptions(), NoSolverStatistics())
