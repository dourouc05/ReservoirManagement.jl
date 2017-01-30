facts("Reservoir data structures: reservoir") do
  n_weeks = 52
  n_scenarios = 5
  year = 1999
  envflow = .02m^3/s
  maxflow = 5.m^3/s

  name1 = "RivN1"
  name2 = "RivN2"
  name3 = "RivD3"
  name4 = "RivD4"

  nr1 = gen_river(name1, n_weeks, year, n_scenarios)
  nr2 = gen_river(name2, n_weeks, year, n_scenarios)
  dr1 = gen_river(name3, n_weeks, year, n_scenarios, envflow, maxflow)
  dr2 = gen_river(name4, n_weeks, year, n_scenarios, envflow, maxflow)

  nrs = [nr1, nr2]
  nrs_names = [name1, name2]
  drs = [dr1, dr2]
  drs_names = [name3, name4]
  rs = [nr1, nr2, dr1, dr2]
  rs_names = [name1, name2, name3, name4]

  r_name = "Reservoir: as simple as possible"
  r_min_capacity = 2.e8m^3
  r_max_capacity = 8.e8m^3

  const a = 1.
  const a_discharge = a * 1.m^3/s
  const a_power = a * 1.W
  const a_height = a * 1.m
  const b = 2.
  const b_discharge = b * 1.m^3/s
  const b_power = b * 1.W
  const c = 3.

  function reservoirTestShortVersions(r::Reservoir)
    outputIds = [:spillway, :penstock, :bottomOutlets]
    purposeIds = [:drinkingWater, :environmentalFlow]
    purposeTypes = [WaterWithdrawalPurpose, DeterministicWaterWithdrawalPurpose, DeterministicDrinkingWater,
                    DeterministicEnvironmentalFlow, DeterministicHydropower]

    @fact name(r) --> getName(r)
    @fact variant(r) --> getVariant(r)
    @fact minCapacity(r) --> getMinimumCapacity(r)
    @fact maxCapacity(r) --> getMaximumCapacity(r)
    @fact penstockHydropower(r) --> getPenstockHydropower(r)
    @fact hydropowerUnits(r) --> getHydropowerUnits(r)
    @fact outputs(r) --> getOutputs(r)
    @fact outputsId(r) --> getOutputIds(r)
    @fact bathymetry(r) --> getBathymetry(r)

    # Getter with an ID: return a Nullable (the other option: throw an error when there is nothing; hard to
    # ensure as a generic test).
    for s in outputIds; @fact output(r, s, nullable=true) --> exactly(getOutput(r, s, nullable=true)); end
    for s in purposeIds; @fact purpose(r, s, nullable=true) --> exactly(getPurpose(r, s, nullable=true)); end

    # Getter with a type: returns arrays, hence no exactly().
    for pt in purposeTypes; @fact purpose(r, pt) --> getPurpose(r, pt); end
  end

  function reservoirTestConsistency(r::Reservoir)
    @fact ReservoirManagement.countScenariosWithConsistency(r) --> countScenarios(r)
    @fact ReservoirManagement.countTimeStepsWithConsistency(r) --> countTimeSteps(r)
    @fact ReservoirManagement.getPeriodWithConsistency(r) --> getPeriod(r)
  end

  function reservoirTestRivers(r::Reservoir, all_rivers, all_natural_rivers, all_diverted_rivers,
                               all_rivers_names, all_natural_rivers_names, all_diverted_rivers_names)
    # Iterators.
    coll_all_rivers = collect(rivers(r))
    coll_all_natural_rivers = collect(naturalRivers(r))
    coll_all_diverted_rivers = collect(divertedRivers(r))
    for riv in all_rivers;          @fact in(riv, coll_all_rivers)          --> true; end
    for riv in all_natural_rivers;  @fact in(riv, coll_all_natural_rivers)  --> true; end
    for riv in all_diverted_rivers; @fact in(riv, coll_all_diverted_rivers) --> true; end

    # Nice iterators.
    each_all_rivers = collect(eachriver(r))
    each_all_natural_rivers = collect(eachnaturalriver(r))
    each_all_diverted_rivers = collect(eachdivertedriver(r))
    for riv in all_rivers;          @fact in(riv, coll_all_rivers)          --> true; end
    for riv in all_natural_rivers;  @fact in(riv, coll_all_natural_rivers)  --> true; end
    for riv in all_diverted_rivers; @fact in(riv, coll_all_diverted_rivers) --> true; end

    # Collectives.
    for riv in all_rivers;          @fact in(riv, getRivers(r))         --> true; end
    for riv in all_natural_rivers;  @fact in(riv, getNaturalRivers(r))  --> true; end
    for riv in all_diverted_rivers; @fact in(riv, getDivertedRivers(r)) --> true; end

    # Raw indexing.
    for i in 1:length(all_natural_rivers);  @fact in(getNaturalRiver(r, i),  all_natural_rivers)  --> true; end
    for i in 1:length(all_diverted_rivers); @fact in(getDivertedRiver(r, i), all_diverted_rivers) --> true; end

    # Name indexing.
    @fact getRiver(r, "Nice name for a river") --> isnull
    @fact hasRiver(r, "Nice name for a river") --> false
    for n in all_rivers_names; @fact getRiver(r, n) --> not(isnull); end
    for n in all_rivers_names; @fact hasRiver(r, n) --> true; end
    for n in all_rivers_names; @fact in(get(getRiver(r, n)), all_rivers) --> true; end

    # Name retrieving.
    for n in all_rivers_names;          @fact in(n, getRiverNames(r))         --> true; end
    for n in all_natural_rivers_names;  @fact in(n, getNaturalRiverNames(r))  --> true; end
    for n in all_diverted_rivers_names; @fact in(n, getDivertedRiverNames(r)) --> true; end
  end

  function reservoirTestScenarioIterator(r::Reservoir)
    for s in eachscenario(r)
      @fact countScenarios(s) --> 1
      for riv in rivers(s)
        actual = getScenario(riv, 1)
        expecteds = [getScenario(get(getRiver(r, getName(riv))), sp) for sp in 1:countScenarios(r)]
        @fact in(actual, expecteds) --> true
      end
    end
  end

  context("Error cases") do
    # No input: neither rivers nor extraction.
    @fact_throws Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity))

    # Duplicate IDs.
    p = DeterministicDrinkingWater(a_discharge)
    @fact_throws Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), purposes=[p, p])
  end

  context("Simplest case") do
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[nr1, nr2], rivers_diverted=[dr1, dr2])

    # Basic getters.
    @fact getName(r) --> r_name
    @fact getMinimumCapacity(r) --> ReservoirManagement._from_unitful(r_min_capacity)
    @fact getMaximumCapacity(r) --> ReservoirManagement._from_unitful(r_max_capacity)
    @fact getPenstockHydropower(r) --> 0.
    @fact getHydropowerUnits(r) --> []
    @fact getPurposes(r) --> []
    @fact getPurposeIds(r) --> []
    @fact getOutputs(r) --> []
    @fact getOutputIds(r) --> []
    @fact getBathymetry(r) --> isempty

    reservoirTestShortVersions(r)
    reservoirTestConsistency(r)
    reservoirTestRivers(r, rs, nrs, drs, rs_names, nrs_names, drs_names)
    reservoirTestScenarioIterator(r)

    # Find purposes and outputs.
    @fact_throws hasPurpose(r, Reservoir)
    @fact hasPurpose(r, WaterWithdrawalPurpose) --> false
    @fact hasPurpose(r, DeterministicWaterWithdrawalPurpose) --> false
    @fact hasPurpose(r, DeterministicDrinkingWater) --> false
    @fact hasPurpose(r, DeterministicEnvironmentalFlow) --> false
    @fact hasPurpose(r, DeterministicHydropower) --> false

    @fact_throws getPurpose(r, Reservoir)
    @fact getPurpose(r, WaterWithdrawalPurpose) --> isempty
    @fact getPurpose(r, DeterministicWaterWithdrawalPurpose) --> isempty
    @fact getPurpose(r, DeterministicDrinkingWater) --> isempty
    @fact getPurpose(r, DeterministicEnvironmentalFlow) --> isempty
    @fact getPurpose(r, DeterministicHydropower) --> isempty

    @fact hasPurpose(r, :drinkingWater) --> false
    @fact hasPurpose(r, :environmentalFlow) --> false
    @fact_throws getPurpose(r, :drinkingWater, nullable=false)
    @fact_throws getPurpose(r, :environmentalFlow, nullable=false)
    @fact getPurpose(r, :drinkingWater, nullable=true) --> isnull
    @fact getPurpose(r, :environmentalFlow, nullable=true) --> isnull

    @fact hasOutput(r, :spillway) --> false
    @fact hasOutput(r, :penstock) --> false
    @fact hasOutput(r, :bottomOutlets) --> false
    @fact_throws getOutput(r, :spillway, nullable=false)
    @fact_throws getOutput(r, :penstock, nullable=false)
    @fact_throws getOutput(r, :bottomOutlets, nullable=false)
    @fact getOutput(r, :spillway, nullable=true) --> isnull
    @fact getOutput(r, :penstock, nullable=true) --> isnull
    @fact getOutput(r, :bottomOutlets, nullable=true) --> isnull

    @fact copy(r) --> r

    # Counting things.
    @fact countScenarios(r) --> n_scenarios
    @fact countTimeSteps(r) --> n_weeks
    @fact getPeriod(r) --> Week(1)
  end

  context("Simplest case, agregating") do
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[nr1, nr2], rivers_diverted=[dr1, dr2])

    @fact_throws totalInflow(r) # Multiple scenarios.
    for s in 1:n_scenarios
      # Inflow as arrays.
      actual = totalInflow(r, s, asTimeSeries=false)
      expected = sum([getScenario(riv, s) for riv in rs])
        @fact actual --> expected

      # Inflow as time series.
      actual = totalInflow(r, s, asTimeSeries=true)
      nr1_sc = getScenarios(nr1)[1]
      expected = TimeArray(timestamp(nr1_sc), sum([getScenario(riv, s) for riv in rs]), colnames(nr1_sc), meta(nr1_sc))
      @fact timestamp(actual) --> timestamp(expected)
      @fact values(actual) --> values(expected)
      @fact colnames(actual) --> colnames(expected)
      @fact meta(actual) --> meta(expected)
    end
  end

  context("Simplest case, altering: deleting a scenario") do
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[nr1, nr2], rivers_diverted=[dr1, dr2])
    r2 = removeScenarios(r, IntSet(1))

    @fact getName(r2) --> r_name
    @fact getMinimumCapacity(r2) --> ReservoirManagement._from_unitful(r_min_capacity)
    @fact getMaximumCapacity(r2) --> ReservoirManagement._from_unitful(r_max_capacity)
    @fact getPenstockHydropower(r2) --> 0.
    @fact getHydropowerUnits(r2) --> []
    @fact getPurposes(r2) --> []
    @fact getPurposeIds(r2) --> []
    @fact getOutputs(r2) --> []
    @fact getOutputIds(r2) --> []
    @fact getBathymetry(r2) --> isempty

    reservoirTestShortVersions(r2)
    reservoirTestConsistency(r2)
    # Cannot apply reservoirTestRivers: the rivers have changed, they have fewer scenarios!
    reservoirTestScenarioIterator(r)

    @fact countScenarios(r2) --> n_scenarios - 1
    @fact countTimeSteps(r2) --> n_weeks
    @fact getPeriod(r2) --> Week(1)

    for rn in rs_names
      for s in 1:n_scenarios - 1
        @fact getScenario(get(getRiver(r2, rn)), s) --> getScenario(get(getRiver(r, rn)), s + 1)
      end
    end

    # Other APIs: don't directly give an IntSet.
    r3 = removeScenarios(r, 1)
    r4 = removeScenarios(r, [1])
    for rn in rs_names
      for s in 1:n_scenarios - 1
        @fact getScenario(get(getRiver(r3, rn)), s) --> getScenario(get(getRiver(r, rn)), s + 1)
        @fact getScenario(get(getRiver(r4, rn)), s) --> getScenario(get(getRiver(r, rn)), s + 1)
      end
    end
  end

  context("Simplest case, altering: deleting multiple scenarios") do
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[nr1, nr2], rivers_diverted=[dr1, dr2])
    r2 = removeScenarios(r, IntSet([1, 3, 5]))

    @fact getName(r2) --> r_name
    @fact getMinimumCapacity(r2) --> ReservoirManagement._from_unitful(r_min_capacity)
    @fact getMaximumCapacity(r2) --> ReservoirManagement._from_unitful(r_max_capacity)
    @fact getPenstockHydropower(r2) --> 0.
    @fact getHydropowerUnits(r2) --> []
    @fact getPurposes(r2) --> []
    @fact getPurposeIds(r2) --> []
    @fact getOutputs(r2) --> []
    @fact getOutputIds(r2) --> []
    @fact getBathymetry(r2) --> isempty

    reservoirTestShortVersions(r2)
    reservoirTestConsistency(r2)
    # Cannot apply reservoirTestRivers: the rivers have changed, they have fewer scenarios!
    reservoirTestScenarioIterator(r)

    @fact countScenarios(r2) --> n_scenarios - 3
    @fact countTimeSteps(r2) --> n_weeks
    @fact getPeriod(r2) --> Week(1)

    for rn in rs_names
      @fact getScenario(get(getRiver(r2, rn)), 1) --> getScenario(get(getRiver(r, rn)), 2)
      @fact getScenario(get(getRiver(r2, rn)), 2) --> getScenario(get(getRiver(r, rn)), 4)
    end

    # Other APIs: don't directly give an IntSet.
    r3 = removeScenarios(r, [1, 3, 5])
    for rn in rs_names
      @fact getScenario(get(getRiver(r3, rn)), 1) --> getScenario(get(getRiver(r, rn)), 2)
      @fact getScenario(get(getRiver(r3, rn)), 2) --> getScenario(get(getRiver(r, rn)), 4)
    end
  end

  context("Simplest case, altering: keeping one time step") do
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[nr1, nr2], rivers_diverted=[dr1, dr2])
    chosen_t_step = 8
    r2 = keepTimeSteps(r, chosen_t_step)

    @fact getName(r2) --> r_name
    @fact getMinimumCapacity(r2) --> ReservoirManagement._from_unitful(r_min_capacity)
    @fact getMaximumCapacity(r2) --> ReservoirManagement._from_unitful(r_max_capacity)
    @fact getPenstockHydropower(r2) --> 0.
    @fact getHydropowerUnits(r2) --> []
    @fact getPurposes(r2) --> []
    @fact getPurposeIds(r2) --> []
    @fact getOutputs(r2) --> []
    @fact getOutputIds(r2) --> []
    @fact getBathymetry(r2) --> isempty

    reservoirTestShortVersions(r2)
    reservoirTestConsistency(r2)
    # Cannot apply reservoirTestRivers: the rivers have changed, they have fewer time steps!
    reservoirTestScenarioIterator(r)

    @fact countScenarios(r2) --> n_scenarios
    @fact countTimeSteps(r2) --> 1
    @fact getPeriod(r2) --> Week(1)

    for rn in rs_names
      for s in 1:n_scenarios
        @fact getScenario(get(getRiver(r2, rn)), s, 1) --> getScenario(get(getRiver(r, rn)), s, chosen_t_step)
      end
    end
  end

  context("Simplest case, altering: keeping multiple time steps") do
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity), rivers_in=[nr1, nr2], rivers_diverted=[dr1, dr2])
    first_chosen_t_step = 8
    last_chosen_t_step = 42
    n_chosen_t_steps = length(collect(first_chosen_t_step:last_chosen_t_step))
    r2 = keepTimeSteps(r, first_chosen_t_step:last_chosen_t_step)

    @fact getName(r2) --> r_name
    @fact getMinimumCapacity(r2) --> ReservoirManagement._from_unitful(r_min_capacity)
    @fact getMaximumCapacity(r2) --> ReservoirManagement._from_unitful(r_max_capacity)
    @fact getPenstockHydropower(r2) --> 0.
    @fact getHydropowerUnits(r2) --> []
    @fact getPurposes(r2) --> []
    @fact getPurposeIds(r2) --> []
    @fact getOutputs(r2) --> []
    @fact getOutputIds(r2) --> []
    @fact getBathymetry(r2) --> isempty

    reservoirTestShortVersions(r2)
    reservoirTestConsistency(r2)
    # Cannot apply reservoirTestRivers: the rivers have changed, they have fewer time steps!
    reservoirTestScenarioIterator(r)

    @fact countScenarios(r2) --> n_scenarios
    @fact countTimeSteps(r2) --> n_chosen_t_steps
    @fact getPeriod(r2) --> Week(1)

    for rn in rs_names
      for s in 1:n_scenarios
        for t in 1:n_chosen_t_steps
          @fact getScenario(get(getRiver(r2, rn)), s, t) --> getScenario(get(getRiver(r, rn)), s, first_chosen_t_step + t - 1)
        end
      end
    end
  end

  context("Only hydropower") do
    hp = HydropowerUnit(max_discharge=a_discharge, max_power=b_power, efficiency=1./c)
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity),
                  rivers_in=[nr1, nr2], rivers_diverted=[dr1, dr2],
                  hydropower=[hp])

    # Basic getters.
    @fact getName(r) --> r_name
    @fact getMinimumCapacity(r) --> ReservoirManagement._from_unitful(r_min_capacity)
    @fact getMaximumCapacity(r) --> ReservoirManagement._from_unitful(r_max_capacity)
    @fact getPenstockHydropower(r) --> ReservoirManagement._from_unitful(a_discharge)
    @fact getHydropowerUnits(r) --> [hp]
    @fact getPurposes(r) --> [DeterministicHydropower(hp)]
    @fact getPurposeIds(r) --> [getPurposeId(hp)]
    @fact getOutputs(r) --> [HydropowerDamOutput(hp)]
    @fact getOutputIds(r) --> [getDamOutputId(hp)]
    @fact getBathymetry(r) --> isempty

    reservoirTestShortVersions(r)
    reservoirTestConsistency(r)
    reservoirTestRivers(r, rs, nrs, drs, rs_names, nrs_names, drs_names)
    reservoirTestScenarioIterator(r)

    # Find purposes and outputs.
    @fact_throws hasPurpose(r, Reservoir)
    @fact hasPurpose(r, WaterWithdrawalPurpose) --> true
    @fact hasPurpose(r, DeterministicWaterWithdrawalPurpose) --> true
    @fact hasPurpose(r, DeterministicDrinkingWater) --> false
    @fact hasPurpose(r, DeterministicEnvironmentalFlow) --> false
    @fact hasPurpose(r, DeterministicHydropower) --> true

    @fact_throws getPurpose(r, Reservoir)
    @fact getPurpose(r, WaterWithdrawalPurpose) --> [DeterministicHydropower(hp)]
    @fact getPurpose(r, DeterministicWaterWithdrawalPurpose) --> [DeterministicHydropower(hp)]
    @fact getPurpose(r, DeterministicDrinkingWater) --> isempty
    @fact getPurpose(r, DeterministicEnvironmentalFlow) --> isempty
    @fact getPurpose(r, DeterministicHydropower) --> [DeterministicHydropower(hp)]

    @fact hasPurpose(r, :drinkingWater) --> false
    @fact hasPurpose(r, :environmentalFlow) --> false
    @fact_throws getPurpose(r, :drinkingWater, nullable=false)
    @fact_throws getPurpose(r, :environmentalFlow, nullable=false)
    @fact getPurpose(r, :drinkingWater, nullable=true) --> isnull
    @fact getPurpose(r, :environmentalFlow, nullable=true) --> isnull

    @fact hasOutput(r, :spillway) --> false
    @fact hasOutput(r, :penstock) --> false
    @fact hasOutput(r, :bottomOutlets) --> false
    @fact_throws getOutput(r, :spillway, nullable=false)
    @fact_throws getOutput(r, :penstock, nullable=false)
    @fact_throws getOutput(r, :bottomOutlets, nullable=false)
    @fact getOutput(r, :spillway, nullable=true) --> isnull
    @fact getOutput(r, :penstock, nullable=true) --> isnull
    @fact getOutput(r, :bottomOutlets, nullable=true) --> isnull

    @fact hasOutput(r, getDamOutputId(hp)) --> true
    @fact getOutput(r, getDamOutputId(hp)) --> HydropowerDamOutput(hp)

    @fact copy(r) --> r

    # Counting things.
    @fact countScenarios(r) --> n_scenarios
    @fact countTimeSteps(r) --> n_weeks
    @fact getPeriod(r) --> Week(1)
  end

  context("Purposes and outputs") do
    p_dw = DeterministicDrinkingWater(a_discharge)
    p_ef = DeterministicEnvironmentalFlow(b_discharge)
    ps = Purpose[p_dw, p_ef]
    o_ep = ConstantDamOutput(:bottomOutlets, a_discharge)
    o_sw = ConstantDamOutput(:spillway, b_discharge)
    os = DamOutput[o_ep, o_sw]
    r = Reservoir(name=r_name, capacity=(r_min_capacity, r_max_capacity),
                  rivers_in=[nr1, nr2], rivers_diverted=[dr1, dr2],
                  purposes=ps, outputs=os)

    # Basic getters.
    @fact getName(r) --> r_name
    @fact getMinimumCapacity(r) --> ReservoirManagement._from_unitful(r_min_capacity)
    @fact getMaximumCapacity(r) --> ReservoirManagement._from_unitful(r_max_capacity)
    @fact getPenstockHydropower(r) --> 0.
    @fact getHydropowerUnits(r) --> isempty
    @fact getPurposes(r) --> ps
    @fact getPurposeIds(r) --> [:drinkingWater, :environmentalFlow]
    @fact getOutputs(r) --> os
    @fact getOutputIds(r) --> [:bottomOutlets, :spillway]
    @fact getBathymetry(r) --> isempty

    reservoirTestShortVersions(r)
    reservoirTestConsistency(r)
    reservoirTestRivers(r, rs, nrs, drs, rs_names, nrs_names, drs_names)
    reservoirTestScenarioIterator(r)

    # Find purposes and outputs.
    @fact_throws hasPurpose(r, Reservoir)
    @fact hasPurpose(r, WaterWithdrawalPurpose) --> true
    @fact hasPurpose(r, DeterministicWaterWithdrawalPurpose) --> true
    @fact hasPurpose(r, DeterministicDrinkingWater) --> true
    @fact hasPurpose(r, DeterministicEnvironmentalFlow) --> true
    @fact hasPurpose(r, DeterministicHydropower) --> false

    @fact_throws getPurpose(r, Reservoir)
    @fact getPurpose(r, WaterWithdrawalPurpose) --> ps
    @fact getPurpose(r, DeterministicWaterWithdrawalPurpose) --> ps
    @fact getPurpose(r, DeterministicDrinkingWater) --> Purpose[p_dw]
    @fact getPurpose(r, DeterministicEnvironmentalFlow) --> Purpose[p_ef]
    @fact getPurpose(r, DeterministicHydropower) --> isempty

    @fact hasPurpose(r, :drinkingWater) --> true
    @fact hasPurpose(r, :environmentalFlow) --> true
    @fact getPurpose(r, :drinkingWater, nullable=false) --> p_dw
    @fact getPurpose(r, :environmentalFlow, nullable=false) --> p_ef
    @fact getPurpose(r, :drinkingWater, nullable=true) --> exactly(Nullable(p_dw))
    @fact getPurpose(r, :environmentalFlow, nullable=true) --> exactly(Nullable(p_ef))

    @fact hasOutput(r, :spillway) --> true
    @fact hasOutput(r, :penstock) --> false
    @fact hasOutput(r, :bottomOutlets) --> true
    @fact getOutput(r, :spillway, nullable=false) --> o_sw
    @fact_throws getOutput(r, :penstock, nullable=false)
    @fact getOutput(r, :bottomOutlets, nullable=false) --> o_ep
    @fact getOutput(r, :spillway, nullable=true) --> exactly(Nullable(o_sw))
    @fact getOutput(r, :penstock, nullable=true) --> isnull
    @fact getOutput(r, :bottomOutlets, nullable=true) --> exactly(Nullable(o_ep))

    @fact copy(r) --> r

    # Counting things.
    @fact countScenarios(r) --> n_scenarios
    @fact countTimeSteps(r) --> n_weeks
    @fact getPeriod(r) --> Week(1)
  end
end
