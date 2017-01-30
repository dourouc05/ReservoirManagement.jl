silentSolver() = GurobiSolver(OutputFlag=0) # TODO: how to silence in a solver-agnostic way? The trick given in JuMP FAQ does not work...

facts("Solvers") do
  const year = 2000
  const a = 1.
  const a_discharge = a * 1.m^3/s
  const b = 2.
  const b_discharge = b * 1.m^3/s

  const riv_name = "Riv"
  const r_name = "Reservoir: tests for solvers"
  const r_min_capacity = 2.e8m^3
  const r_max_capacity = 8.e8m^3

  const time_horizon = 5
  const year = 2000
  const period = Day(1)

  context("Helpers: purpose (both purposes, dependence on period)") do
    purposes = DeterministicPurposes(drinkingWater=a_discharge, environmentalFlow=b_discharge)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Day(1)) --> roughly((a + b) * 86400.e-6)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Day(2)) --> roughly((a + b) * 86400.e-6 * 2)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Week(1)) --> roughly((a + b) * 86400.e-6 * 7)
  end

  context("Helpers: purpose (both purposes, same role for both)") do
    # Both purposes have the same impact: switch the discharges, should get the same results as before.
    purposes = DeterministicPurposes(drinkingWater=b_discharge, environmentalFlow=a_discharge)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Day(1)) --> roughly((a + b) * 86400.e-6)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Day(2)) --> roughly((a + b) * 86400.e-6 * 2)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Week(1)) --> roughly((a + b) * 86400.e-6 * 7)
  end

  context("Helpers: purpose (just drinking water)") do
    purposes = DeterministicPurposes(drinkingWater=a_discharge)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Day(1)) --> roughly(a * 86400.e-6)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Day(2)) --> roughly(a * 86400.e-6 * 2)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Week(1)) --> roughly(a * 86400.e-6 * 7)
  end

  context("Helpers: purpose (just environmental flow)") do
    purposes = DeterministicPurposes(drinkingWater=b_discharge)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Day(1)) --> roughly(b * 86400.e-6)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Day(2)) --> roughly(b * 86400.e-6 * 2)
    @fact hlp_JuMP_purposes_discharge(Model(), purposes, Week(1)) --> roughly(b * 86400.e-6 * 7)
  end

  context("Helpers: output (nonlinear modes)") do
    # Recognised (non)linear modes.
    @fact _hlp_outputs_allowed_nonlinear_mode(:None) --> nothing
    @fact _hlp_outputs_allowed_nonlinear_mode(:LinearApprox) --> nothing
    @fact _hlp_outputs_allowed_nonlinear_mode(:PiecewiseLinear) --> nothing
    @fact _hlp_outputs_allowed_nonlinear_mode(:Cone) --> nothing
    @fact_throws _hlp_outputs_allowed_nonlinear_mode(:SomethingElseThatIsNotRecognised)
  end

  function outputHelpersTestLinear(period::Period, length::Int)
    riv = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(length, year), gen_scenario(length), ["Discharge (m³/s)"], period)])
    outputs = ConstantDamOutputs(bottomOutlets=a_discharge, spillway=b_discharge)
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[riv], outputs=outputs)

    # Smallest optimisation model so that this works.
    m = Model(solver=silentSolver())
    @variable(m, output[1:length] >= 0)
    @variable(m, ReservoirManagement._from_unitful(r_min_capacity) <= volume[1:length] <= ReservoirManagement._from_unitful(r_max_capacity))
    @objective(m, Max, sum(output))
    for t in 1:length
      hlp_JuMP_outputs_constraint(m, output[t], volume[t], r, nonlinear_mode=:None)
    end
    status = solve(m)

    # Check the output values: must be the maximum possible to release (the volume is free, there is zero constraint on it).
    @fact status --> :Optimal
    @fact getvalue(sum(output)) --> roughly((a + b) * 86400.e-6 * length * Dates.days(period))
  end

  context("Helpers: output (linear)") do
    # Test the constraints are the expected ones: maximise the output, that must be the input value. It depends on the
    # period and complete duration of optimisation.
    outputHelpersTestLinear(Day(1), 1)
    outputHelpersTestLinear(Day(2), 7)
    outputHelpersTestLinear(Week(1), 12)
  end

  function outputHelpersTestLinearConditional(period::Period, length::Int)
    # bathymetry: h = 280 [m] + V [10^6 m^3] * .1 [m / 10^6 m^3]
    # Bounds: V in [r_min_capacity, r_max_capacity], [m^3] --> 200 to 800 [10^6 m^3].
    # Hence for reservoir level: [300, 360] m. Condition in the middle
    geom = LinearReservoirBathymetry(280., .1)
    riv = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(length, year), gen_scenario(length), ["Discharge (m³/s)"], period)])

    spillway_height = 330.
    r_cap_spillway = computeVolume(geom, spillway_height)
    outputs = ConstantDamOutputs(bottomOutlets=a_discharge, spillway=b_discharge, spillway_minLevel=spillway_height * 1.m)

    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[riv], outputs=outputs, bathymetry=geom)

    ## Disable the spillway.
    # Smallest optimisation model so that this works: completely disable the spillway by working on the maximum volume.
    m_1 = Model(solver=silentSolver())
    @variable(m_1, output_1[1:length] >= 0)
    @variable(m_1, ReservoirManagement._from_unitful(r_min_capacity) <= volume_1[1:length] <= r_cap_spillway - .1) # Need a small epsilon to be sure to be below the spillway!
    @objective(m_1, Max, sum(output_1))
    for t in 1:length
      hlp_JuMP_outputs_constraint(m_1, output_1[t], volume_1[t], r, nonlinear_mode=:None)
    end
    status_1 = solve(m_1)

    # Check the output values: must be the maximum possible to release (the volume is free,
    # there is zero constraint on it, except to disable the spillway).
    @fact status_1 --> :Optimal
    @fact getvalue(sum(output_1)) --> roughly(a * 86400.e-6 * length * Dates.days(period))

    ## Force enable the spillway.
    m_2 = Model(solver=silentSolver())
    @variable(m_2, output_2[1:length] >= 0)
    @variable(m_2, r_cap_spillway <= volume_2[1:length] <= ReservoirManagement._from_unitful(r_max_capacity)) # With the objective function, the spillway is always active.
    @objective(m_2, Max, sum(output_2))
    for t in 1:length
      hlp_JuMP_outputs_constraint(m_2, output_2[t], volume_2[t], r, nonlinear_mode=:None)
    end
    status_2 = solve(m_2)

    # Check the output values: must be the maximum possible to release (the volume is free,
    # there is zero constraint on it). As a consequence, the volume will be chosen to enable the spillway.
    @fact status_2 --> :Optimal
    @fact getvalue(sum(output_2)) --> roughly((a + b) * 86400.e-6 * length * Dates.days(period))

    ## Spillway free (but release maximised).
    m_3 = Model(solver=silentSolver())
    @variable(m_3, output_3[1:length] >= 0)
    @variable(m_3, ReservoirManagement._from_unitful(r_min_capacity) <= volume_3[1:length] <= ReservoirManagement._from_unitful(r_max_capacity))
    @objective(m_3, Max, sum(output_3))
    for t in 1:length
      hlp_JuMP_outputs_constraint(m_3, output_3[t], volume_3[t], r, nonlinear_mode=:None)
    end
    status_3 = solve(m_3)

    # Check the output values: must be the maximum possible to release (the volume is free,
    # there is zero constraint on it, except to force the spillway to be operational).
    @fact status_3 --> :Optimal
    @fact getvalue(sum(output_3)) --> roughly((a + b) * 86400.e-6 * length * Dates.days(period))
  end

  context("Helpers: output (linear, conditional [reservoir level])") do
    # Test the constraints are the expected ones: maximise the output, that must be the input value. It depends on the
    # period and complete duration of optimisation.
    outputHelpersTestLinearConditional(Day(1), 1)
    outputHelpersTestLinearConditional(Day(2), 7)
    outputHelpersTestLinearConditional(Week(1), 12)

    # If the maximum reservoir level is higher than the condition, show a warning.
    # TODO: no test framework allows this for now, it seems! 
  end

  context("Helpers: inflow (one natural river)") do
    # Simplest case: one natural river.
    riv1 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[riv1])

    for t in 1:time_horizon
      @fact hlp_JuMP_inflow_discharge(Model(), r, 1, t) --> roughly(getScenario(riv1, 1, t) * Dates.days(period))
    end
  end

  context("Helpers: inflow (many natural rivers)") do
    riv1 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    riv2 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    riv3 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[riv1, riv2, riv3])

    for t in 1:time_horizon
      @fact hlp_JuMP_inflow_discharge(Model(), r, 1, t) --> roughly((getScenario(riv1, 1, t) + getScenario(riv2, 1, t) + getScenario(riv3, 1, t)) * Dates.days(period))
    end
  end

  context("Helpers: inflow (diverted rivers, imposed)") do
    riv1 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    riv2 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    riv3 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    driv1 = DivertedRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)], environmental_flow=a_discharge, maximum_flow=b_discharge)
    driv2 = DivertedRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)], environmental_flow=a_discharge, maximum_flow=b_discharge)
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[riv1, riv2, riv3], rivers_diverted=[driv1, driv2])

    for t in 1:time_horizon
      @fact hlp_JuMP_inflow_discharge(Model(), r, 1, t, imposeDiverted=true) --> roughly((maxAllowableFlow(driv1, 1, t) + maxAllowableFlow(driv2, 1, t) + getScenario(riv1, 1, t) + getScenario(riv2, 1, t) + getScenario(riv3, 1, t)) * Dates.days(period))
    end
  end

  context("Helpers: inflow (diverted rivers, degree of freedom)") do
    # Must optimise, as the inflow becomes an optimisation variable!
    riv1 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    riv2 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    riv3 = NaturalRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)])
    driv1 = DivertedRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)], environmental_flow=a_discharge, maximum_flow=b_discharge)
    driv2 = DivertedRiver(name=riv_name, scenarios=[TimeArray(gen_dates(time_horizon, year), gen_scenario(time_horizon), ["Discharge (m³/s)"], period)], environmental_flow=a_discharge, maximum_flow=b_discharge)
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[riv1, riv2, riv3], rivers_diverted=[driv1, driv2])

    m_m = Model(solver=silentSolver())
    @variable(m_m, v_divert[1:1, 1:time_horizon, getDivertedRiverNames(r)] >= 0)
    @variable(m_m, obj >= 0)
    @objective(m_m, Max, obj)
    @constraint(m_m, obj == sum{hlp_JuMP_inflow_discharge(m_m, r, driv1, 1, t, imposeDiverted=v_divert) + hlp_JuMP_inflow_discharge(m_m, r, driv2, 1, t, imposeDiverted=v_divert), t in 1:time_horizon})
    solve(m_m)

    @fact getvalue(obj) --> roughly(sum([(maxAllowableFlow(driv1, 1, t) + maxAllowableFlow(driv2, 1, t)) * Dates.days(period) for t in 1:time_horizon]))
  end

  context("Evaluation: purpose shortage (feasible solution)") do
    # 2 m^3/s available, withdraw 1 m^3/s, no output limit.
    ts = TimeArray(gen_dates(time_horizon, year), [2.e-6 * 86400 * Dates.days(period) for t in 1:time_horizon], ["Discharge (m³/s)"], period)
    riv = NaturalRiver(name=riv_name, scenarios=[ts])
    purpose = DeterministicPurposes(drinkingWater=1.0m^3/s)
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[riv], purposes=purpose)

    demands = [sum([getNeed(p, period) for p in getPurposes(r)]) for t in 1:time_horizon]
    sol = minimumRuleCurve(r, verbosity=-1)
    eval = purposeShortage(sol, values(ts), verbosity=-1)

    @fact typeof(eval) --> ReservoirSolution
    @fact typeof(eval.solution) --> ShortageSolution
    @fact typeof(eval.solverOptions) --> ShortageOptions
    @fact typeof(eval.solverStatistics) --> ShortageSolverStatistics
    @fact eval.solver --> :ShortageSolver

    @fact eval.solution.infeasibilityCause --> :None
    @fact eval.solution.delta --> roughly(zeros(time_horizon)) # Structural property and actual expected value.
    @fact sum(eval.solution.delta .<= demands) --> length(eval.solution.delta) # Cannot lack more water than is needed, for each time step.

    @fact computePurposeShortageIndex(eval, purpose) --> roughly(0.) # Structural property and actual expected value.
    @fact computePurposeVulnerabilityIndex(eval, purpose) --> roughly(0.) # Structural property and actual expected value.
  end

  context("Evaluation: purpose shortage (shortage)") do
    # 0 m^3/s available, withdraw 2 m^3/s, no output limit, initial condition forced to the minimum reservoir level.
    ts = TimeArray(gen_dates(time_horizon, year), [0. for t in 1:time_horizon], ["Discharge (m³/s)"], period) # 1.e-6 * 86400 * Dates.days(period)
    riv = NaturalRiver(name=riv_name, scenarios=[ts])
    purpose = DeterministicPurposes(drinkingWater=2.m^3/s)
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[riv], purposes=purpose)

    # Compute a rule curve and evaluate it.
    sol = minimumRuleCurve(r, verbosity=-1)
    eval = purposeShortage(sol, values(ts), forceInitialCondition=true, verbosity=-1)

    # Compute data required for the tests.
    minLevel = getSolution(sol)
    rule_curve_delta = [minLevel[t] - (t == 1 ? getMinimumCapacity(r) : minLevel[t - 1]) for t in 1:time_horizon]
    demands = [sum([getNeed(p, period) for p in getPurposes(r)]) + (rule_curve_delta[t] >= 0. ? rule_curve_delta[t] : 0.) for t in 1:time_horizon]

    # Actual tests.
    @fact typeof(eval) --> ReservoirSolution
    @fact typeof(eval.solution) --> ShortageSolution
    @fact typeof(eval.solverOptions) --> ShortageOptions
    @fact typeof(eval.solverStatistics) --> ShortageSolverStatistics
    @fact eval.solver --> :ShortageSolver

    @fact eval.solution.infeasibilityCause --> :Shortage
    @fact sum(eval.solution.delta) --> greater_than(0.1) # Structural property.
    @fact sum(eval.solution.delta) --> roughly(.864) # Actual expected value.
    @fact sum(eval.solution.delta .<= demands) --> length(eval.solution.delta) # Cannot lack more water than is needed, for each time step.

    @fact computePurposeShortageIndex(eval, purpose) --> greater_than(0.1) # Structural property.
    @fact computePurposeShortageIndex(eval, purpose) --> roughly(1.) # Actual expected value.
    @fact computePurposeVulnerabilityIndex(eval, purpose) --> greater_than(0.1) # Structural property.
    @fact computePurposeVulnerabilityIndex(eval, purpose) --> roughly(5.) # Actual expected value.
    # Why does those two indices give very different results? There is no water in, so the lack of water is very large.
    # The first one considers averages over years (for the only year in the problem, 100% of the water is missing).
    # The second one considers every time step (here, five). As the sum happens over the *squares* of ratios, this gives
    # an impressive difference.
  end

  context("Evaluation: purpose shortage (abundance)") do
    # Time series for deriving the rule curve (ts_lo), ensures there is a feasible solution.
    # Check for infeasibility against more extreme scenarios (ts_hi).
    ts_lo = TimeArray(gen_dates(time_horizon, year), [50.e-6 * 86400 * Dates.days(period) for t in 1:time_horizon], ["Discharge (m³/s)"], period) # 10^6 m^3 / period
    ts_hi = TimeArray(gen_dates(time_horizon, year), [50.e-2 * 86400 * Dates.days(period) for t in 1:time_horizon], ["Discharge (m³/s)"], period) # 10^6 m^3 / period
    riv = NaturalRiver(name=riv_name, scenarios=[ts_lo])
    purpose = DeterministicPurposes(drinkingWater=0.01m^3/s)
    outs = ConstantDamOutputs(bottomOutlets=1.m^3/s)
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[riv], purposes=purpose, outputs=outs)

    demands = [sum([getNeed(p, period) for p in getPurposes(r)]) for t in 1:time_horizon]
    sol = minimumRuleCurve(r, useOutputs=getOutputIds(r), verbosity=-1)
    eval = purposeShortage(sol, values(ts_hi), useOutputs=getOutputIds(r), forceInitialCondition=true, verbosity=-1)

    @fact typeof(eval) --> ReservoirSolution
    @fact typeof(eval.solution) --> ShortageSolution
    @fact typeof(eval.solverOptions) --> ShortageOptions
    @fact typeof(eval.solverStatistics) --> ShortageSolverStatistics
    @fact eval.solver --> :ShortageSolver

    @fact eval.solution.infeasibilityCause --> :Surplus
    @fact sum(eval.solution.delta) --> less_than(0.0) # Structural property: surplus, hence negative.
    @fact abs(sum(eval.solution.delta)) --> greater_than(0.1) # Structural property: far from 0.
    @fact abs(sum(eval.solution.delta)) --> roughly(215999.56367) # Actual expected value.
    # No test between eval.solution.delta and demands, as purposes are fulfilled (surplus of water, not lack).

    @fact computePurposeShortageIndex(eval, purpose) --> roughly(0.) # Structural property and actual expected value.
    @fact computePurposeVulnerabilityIndex(eval, purpose) --> roughly(0.) # Structural property and actual expected value.
  end
end
