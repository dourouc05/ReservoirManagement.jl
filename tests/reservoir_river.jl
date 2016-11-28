facts("Reservoir data structures: rivers") do
  const name1 = "River_name 1"
  const name2 = "River_name*2"
  const a = 1.
  const a_discharge = a * 1.m^3/s
  const b = 4. # Value chosen for DivertedRiver, the function getMaximumAllowableFlow
  const b_discharge = b * 1.m^3/s
  const sc1_begin_year = 1999
  const sc1_n_weeks = 5
  const sc1_n_scenarios = 25
  const sc1 = gen_tses(sc1_n_weeks, sc1_begin_year, sc1_n_scenarios)
  const sc2_begin_year = 1804
  const sc2_n_weeks = 5
  const sc2_n_scenarios = 25
  const sc2 = gen_tses(sc2_n_weeks, sc2_begin_year, sc2_n_scenarios)
  const sc_wrong = [TimeArray(timestamp(sc1[1]), values(sc1[1]), colnames(sc1[1]), week)]
  const sc_nan = [TimeArray(timestamp(sc1[1]), values(sc1[1]) * NaN, colnames(sc1[1]), week)]
  const sc_empty = similar(sc1, 0)

  function naturalRiverTestShortVersions(r::NaturalRiver)
    @fact name(r) --> getName(r)
    @fact scenarios(r) --> getScenarios(r)
    @fact scenario(r, 1) --> getScenario(r, 1)
    @fact scenario(r, 1, 1) --> getScenario(r, 1, 1)
    @fact period(r) --> getPeriod(r)
    @fact countScenarios(r) --> countScenariosWithConsistency(r)
    @fact countTimeSteps(r) --> countTimeStepsWithConsistency(r)
    @fact getPeriod(r) --> getPeriodWithConsistency(r)
    @fact typeof(getPeriod(r)) <: Base.Dates.Period --> true
  end

  context("NaturalRiver") do
    r = NaturalRiver(name=name1, scenarios=sc1)
    @fact isempty(r) --> false
    @fact getName(r) --> name1
    @fact getScenarios(r) --> sc1
    @fact countScenarios(r) --> sc1_n_scenarios
    @fact countTimeSteps(r) --> sc1_n_weeks
    @fact getPeriod(r) --> Week(1)
    @fact getScenario(r, 1) --> values(sc1[1])
    @fact getScenario(r, 1, 1) --> values(sc1[1])[1]
    naturalRiverTestShortVersions(r)

    r2 = copy(r, name=name2)
    @fact isempty(r2) --> false
    @fact getName(r2) --> name2
    @fact getScenarios(r2) --> sc1
    @fact countScenarios(r2) --> sc1_n_scenarios
    @fact countTimeSteps(r2) --> sc1_n_weeks
    @fact getPeriod(r2) --> Week(1)
    @fact getScenario(r2, 1) --> values(sc1[1])
    @fact getScenario(r2, 1, 1) --> values(sc1[1])[1]
    naturalRiverTestShortVersions(r2)

    r3 = copy(r2, scenarios=sc2)
    @fact isempty(r3) --> false
    @fact getName(r3) --> name2
    @fact getScenarios(r3) --> sc2
    @fact countScenarios(r3) --> sc2_n_scenarios
    @fact countTimeSteps(r3) --> sc2_n_weeks
    @fact getPeriod(r3) --> Week(1)
    @fact getScenario(r3, 1) --> values(sc2[1])
    @fact getScenario(r3, 1, 1) --> values(sc2[1])[1]
    naturalRiverTestShortVersions(r3)

    r4 = copy(r, name=name2, scenarios=sc2)
    @fact r3 --> exactly(r4)

    @fact_throws NaturalRiver(name=name1, scenarios=sc_wrong)
    @fact_throws NaturalRiver(name=name1, scenarios=sc_nan)
    @fact_throws NaturalRiver(name=name1, scenarios=sc_empty)
  end

  function divertedRiverTestShortVersions(r::DivertedRiver)
    @fact name(r) --> getName(r)
    @fact scenarios(r) --> getScenarios(r)
    @fact scenario(r, 1) --> getScenario(r, 1)
    @fact scenario(r, 1, 1) --> getScenario(r, 1, 1)
    @fact period(r) --> getPeriod(r)
    @fact countScenarios(r) --> countScenariosWithConsistency(r)
    @fact countTimeSteps(r) --> countTimeStepsWithConsistency(r)
    @fact getPeriod(r) --> getPeriodWithConsistency(r)
    @fact typeof(getPeriod(r)) <: Base.Dates.Period --> true
    @fact environmentalFlow(r) --> getEnvironmentalFlow(r)
    @fact maxFlow(r) --> getMaximumFlow(r)
    @fact environmentalFlow(r, Week(1)) --> getEnvironmentalFlow(r, Week(1))
    @fact maxFlow(r, Week(1)) --> getMaximumFlow(r, Week(1))
    @fact maxAllowableFlow(r, 1, 1) --> getMaximumAllowableFlow(r, 1, 1)
  end

  context("DivertedRiver") do
    r = DivertedRiver(name=name1, scenarios=sc1, environmental_flow=a_discharge, maximum_flow=b_discharge)
    @fact isempty(r) --> false
    @fact getName(r) --> name1
    @fact getScenarios(r) --> sc1
    @fact countScenarios(r) --> sc1_n_scenarios
    @fact countTimeSteps(r) --> sc1_n_weeks
    @fact getPeriod(r) --> Week(1)
    @fact getScenario(r, 1) --> values(sc1[1])
    @fact getScenario(r, 1, 1) --> values(sc1[1])[1]
    @fact getEnvironmentalFlow(r) --> _from_unitful(a_discharge)
    @fact getMaximumFlow(r) --> _from_unitful(b_discharge)
    @fact getEnvironmentalFlow(r, Week(1)) --> _from_unitful(a_discharge) * 604800. # 10^6 m^3/s -> 10^6 m^3/week
    @fact getMaximumFlow(r, Week(1)) --> _from_unitful(b_discharge) * 604800. # 10^6 m^3/s -> 10^6 m^3/week
    @fact maxAllowableFlow(r, 1, 1) --> roughly(maximum([0., minimum([getMaximumFlow(r, Week(1)), getScenario(r, 1, 1) - getEnvironmentalFlow(r, Week(1))])]))
    divertedRiverTestShortVersions(r)

    r2 = copy(r, name=name2)
    @fact isempty(r2) --> false
    @fact getName(r2) --> name2
    @fact getScenarios(r2) --> sc1
    @fact countScenarios(r2) --> sc1_n_scenarios
    @fact countTimeSteps(r2) --> sc1_n_weeks
    @fact getPeriod(r2) --> Week(1)
    @fact getScenario(r2, 1) --> values(sc1[1])
    @fact getScenario(r2, 1, 1) --> values(sc1[1])[1]
    @fact getEnvironmentalFlow(r2) --> _from_unitful(a_discharge)
    @fact getMaximumFlow(r2) --> _from_unitful(b_discharge)
    @fact maxAllowableFlow(r2, 1, 1) --> roughly(maximum([0., minimum([getMaximumFlow(r2, Week(1)), getScenario(r2, 1, 1) - getEnvironmentalFlow(r2, Week(1))])]))
    divertedRiverTestShortVersions(r2)

    r3 = copy(r2, scenarios=sc2)
    @fact isempty(r3) --> false
    @fact getName(r3) --> name2
    @fact getScenarios(r3) --> sc2
    @fact countScenarios(r3) --> sc2_n_scenarios
    @fact countTimeSteps(r3) --> sc2_n_weeks
    @fact getPeriod(r3) --> Week(1)
    @fact getScenario(r3, 1) --> values(sc2[1])
    @fact getScenario(r3, 1, 1) --> values(sc2[1])[1]
    @fact getEnvironmentalFlow(r3) --> _from_unitful(a_discharge)
    @fact getMaximumFlow(r3) --> _from_unitful(b_discharge)
    @fact maxAllowableFlow(r3, 1, 1) --> roughly(maximum([0., minimum([getMaximumFlow(r3, Week(1)), getScenario(r3, 1, 1) - getEnvironmentalFlow(r3, Week(1))])]))
    divertedRiverTestShortVersions(r3)

    r4 = copy(r3, environmental_flow=b_discharge)
    @fact isempty(r4) --> false
    @fact getName(r4) --> name2
    @fact getScenarios(r4) --> sc2
    @fact countScenarios(r4) --> sc2_n_scenarios
    @fact countTimeSteps(r4) --> sc2_n_weeks
    @fact getPeriod(r4) --> Week(1)
    @fact getScenario(r4, 1) --> values(sc2[1])
    @fact getScenario(r4, 1, 1) --> values(sc2[1])[1]
    @fact getEnvironmentalFlow(r4) --> _from_unitful(b_discharge)
    @fact getMaximumFlow(r4) --> _from_unitful(b_discharge)
    @fact maxAllowableFlow(r4, 1, 1) --> roughly(maximum([0., minimum([getMaximumFlow(r4, Week(1)), getScenario(r4, 1, 1) - getEnvironmentalFlow(r4, Week(1))])]))
    divertedRiverTestShortVersions(r4)

    r5 = copy(r4, maximum_flow=a_discharge)
    @fact isempty(r5) --> false
    @fact getName(r5) --> name2
    @fact getScenarios(r5) --> sc2
    @fact countScenarios(r5) --> sc2_n_scenarios
    @fact countTimeSteps(r5) --> sc2_n_weeks
    @fact getPeriod(r5) --> Week(1)
    @fact getScenario(r5, 1) --> values(sc2[1])
    @fact getScenario(r5, 1, 1) --> values(sc2[1])[1]
    @fact getEnvironmentalFlow(r5) --> _from_unitful(b_discharge)
    @fact getMaximumFlow(r5) --> _from_unitful(a_discharge)
    @fact maxAllowableFlow(r5, 1, 1) --> roughly(maximum([0., minimum([getMaximumFlow(r5, Week(1)), getScenario(r5, 1, 1) - getEnvironmentalFlow(r5, Week(1))])]))
    divertedRiverTestShortVersions(r5)

    r6 = copy(r, name=name2, scenarios=sc2, environmental_flow=b_discharge, maximum_flow=a_discharge)
    @fact r6 --> exactly(r5)

    @fact_throws DivertedRiver(name=name1, scenarios=sc_wrong, environmental_flow=a_discharge, maximum_flow=b_discharge)
    @fact_throws DivertedRiver(name=name1, scenarios=sc_nan, environmental_flow=a_discharge, maximum_flow=b_discharge)
    @fact_throws DivertedRiver(name=name1, scenarios=sc_empty, environmental_flow=a_discharge, maximum_flow=b_discharge)
  end
end
