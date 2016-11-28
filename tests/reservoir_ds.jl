# Type definitions for test "Purpose: general implementation of getId works for user-defined purposes".
immutable Test1Purpose <: Purpose
  id::Symbol
end
immutable Test2Purpose <: Purpose
  not_id::Symbol
end
getId(p::Test2Purpose) = p.not_id

# Type definitions for test "DamOutput: general implementation of getId works for user-defined releases".
immutable Test1DamOutput <: DamOutput
  id::Symbol
end
immutable Test2DamOutput <: DamOutput
  not_id::Symbol
end
getId(p::Test2DamOutput) = p.not_id

# Actual tests.
facts("Reservoir data structures: basics") do
  function bathymetryTestShortVersions(g::ReservoirBathymetry)
    @fact isPolynomial(g) --> ispolynomial(g)
    @fact isLinear(g) --> islinear(g)
    @fact isQuadratic(g) --> isquadratic(g)
    @fact isPiecewiseLinear(g) --> ispiecewiselinear(g)
    @fact isPWL(g) --> ispiecewiselinear(g)
    @fact isSquareRoot(g) --> issquareroot(g)
  end

  function bathymetryTestValueShortVersions(g::ReservoirBathymetry)
    @fact volume(g, 2.0) --> computeVolume(g, 2.0)
    @fact volume(g, 4.0) --> computeVolume(g, 4.0)
    @fact level(g, 2.0) --> computeLevel(g, 2.0)
    @fact level(g, 4.0) --> computeLevel(g, 4.0)
  end

  const an_id = :a
  const a = 1.
  const a_discharge = a * 1.m^3/s
  const a_power = a * 1.W
  const a_height = a * 1.m
  const b = 2.
  const b_discharge = b * 1.m^3/s
  const b_power = b * 1.W
  const c = 3.
  const xs = [1.1, 2.0, 3.0, 4.0] # xs[1] = 1.1 for QuadraticReservoirBathymetry: avoid root 0 (so xs[1] != 1.0).

  context("bathymetry: NullReservoirBathymetry") do
    g = NullReservoirBathymetry()
    @fact isempty(g) --> true
    for x in xs; @fact_throws level(g, x); end
    for x in xs; @fact_throws volume(g, x); end
    @fact iswhitebox(g) --> false
    # What about isblackbox? No value makes sense...
    @fact ispolynomial(g) --> false
    @fact islinear(g) --> false
    @fact isquadratic(g) --> false
    @fact ispiecewiselinear(g) --> false
    @fact issquareroot(g) --> false
    bathymetryTestShortVersions(g)

    @fact_throws computeLevel(g, 2.0)
    @fact_throws level(g, 2.0)
    @fact_throws computeVolume(g, 2.0)
    @fact_throws volume(g, 2.0)
  end

  context("bathymetry: LinearReservoirBathymetry") do
    g = LinearReservoirBathymetry((a, b))
    @fact isempty(g) --> false
    for x in xs; @fact level(g, x) --> roughly(a + b * x); end
    for x in xs; @fact volume(g, x)  --> roughly((x - a) / b); end
    @fact iswhitebox(g) --> true
    @fact isblackbox(g) --> false
    @fact ispolynomial(g) --> true
    @fact islinear(g) --> true
    @fact isquadratic(g) --> false
    @fact ispiecewiselinear(g) --> false
    @fact issquareroot(g) --> false
    bathymetryTestShortVersions(g)
    bathymetryTestValueShortVersions(g)

    @fact coefficient(g, 1) --> b
    @fact coefficient(g, 0) --> a

    # Other constructors.
    g2 = LinearReservoirBathymetry(a, b)
    @fact g2 --> exactly(g)
  end

  context("bathymetry: QuadraticReservoirBathymetry") do
    g = QuadraticReservoirBathymetry((a, b, c))
    # Roots for this polynomial:
    #     {{v -> 1/3 (-1 - Sqrt[-2 + 3 d])}, {v -> 1/3 (-1 + Sqrt[-2 + 3 d])}}
    @fact isempty(g) --> false
    for x in xs; @fact level(g, x) --> roughly(a + b * x + c * x^2); end
    for x in xs; @fact volume(g, x)  --> roughly((- 1 + sqrt(- 2 + 3 * x))/ 3); end
    @fact iswhitebox(g) --> true
    @fact isblackbox(g) --> false
    @fact ispolynomial(g) --> true
    @fact islinear(g) --> false
    @fact isquadratic(g) --> true
    @fact ispiecewiselinear(g) --> false
    @fact issquareroot(g) --> false
    bathymetryTestShortVersions(g)
    bathymetryTestValueShortVersions(g)

    @fact coefficient(g, 2) --> c
    @fact coefficient(g, 1) --> b
    @fact coefficient(g, 0) --> a

    @fact_throws volume(QuadraticReservoirBathymetry((-1., -2., -3.)), 1.1) # No real root
    @fact_throws volume(QuadraticReservoirBathymetry((1.21, -2.2, 1.)), 1.1) # Two roots.
    @fact_throws volume(QuadraticReservoirBathymetry((1., 2., 3.)), 1.) # Root is zero.

    # Other constructors.
    g2 = QuadraticReservoirBathymetry(a, b, c)
    @fact g2 --> exactly(g)
  end

  context("bathymetry: PiecewiseLinearReservoirBathymetry") do
    @fact_throws PiecewiseLinearReservoirBathymetry([a, b, c], xs[1:3]) # Too many values
    @fact_throws PiecewiseLinearReservoirBathymetry([a, b, c], xs[1:1]) # Too few values

    g = PiecewiseLinearReservoirBathymetry([a, b, c], xs[1:2])
    @fact isempty(g) --> false
    # for x in xs; @fact level(g, x) --> roughly(...); end
    # for x in xs; @fact volume(g, x)  --> roughly(...); end
    @fact iswhitebox(g) --> true
    @fact isblackbox(g) --> false
    @fact ispolynomial(g) --> false
    @fact islinear(g) --> false
    @fact isquadratic(g) --> false
    @fact ispiecewiselinear(g) --> true
    @fact issquareroot(g) --> false
    bathymetryTestShortVersions(g)
    # TODO: bathymetryTestValueShortVersions(g)

    # TODO: values.
  end

  context("bathymetry: SquareRootReservoirBathymetry") do
    g = SquareRootReservoirBathymetry(a, b)
    @fact isempty(g) --> false
    for x in xs; @fact level(g, x) --> roughly(a + sqrt(b * x)); end
    for x in xs; @fact volume(g, x)  --> roughly((x - a)^2 / b); end
    @fact iswhitebox(g) --> true
    @fact isblackbox(g) --> false
    @fact ispolynomial(g) --> false
    @fact islinear(g) --> false
    @fact isquadratic(g) --> false
    @fact ispiecewiselinear(g) --> false
    @fact issquareroot(g) --> true
    bathymetryTestShortVersions(g)
    bathymetryTestValueShortVersions(g)

    # Other constructors.
    g2 = SquareRootReservoirBathymetry((a, b))
    @fact g2 --> exactly(g)
  end

  context("HydropowerUnit") do
    @fact_throws HydropowerUnit(max_discharge=a_discharge, max_power=b_power, efficiency=2.) # Wrong efficiency

    hp = HydropowerUnit(max_discharge=a_discharge, max_power=b_power, efficiency=1./c)
    @fact getId(hp) --> not(getPurposeId(hp))
    @fact getId(hp) --> not(getDamOutputId(hp))
    @fact getPurposeId(hp) --> not(getDamOutputId(hp))
    @fact getMaximumDischarge(hp) --> _from_unitful(a_discharge)
    @fact getMaximumPower(hp) --> _from_unitful(b_power)
    @fact getEfficiency(hp) --> 1./c

    @fact id(hp) --> getId(hp)
    @fact purposeId(hp) --> getPurposeId(hp)
    @fact damOutputId(hp) --> getDamOutputId(hp)
    @fact maxDischarge(hp) --> getMaximumDischarge(hp)
    @fact maxPower(hp) --> getMaximumPower(hp)
    @fact efficiency(hp) --> getEfficiency(hp)
  end

  context("Purpose: general implementation of methods for user-defined purposes") do
    # Use the existing method, based on the id field.
    t = Test1Purpose(an_id)
    @fact getId(t) --> an_id
    @fact id(t) --> an_id
    @fact isWaterWithdrawal(t) --> false
    @fact isDeterministic(t) --> false

    # Use the user-defined getId.
    t = Test2Purpose(an_id)
    @fact getId(t) --> an_id
    @fact id(t) --> an_id
    @fact isWaterWithdrawal(t) --> false
    @fact isDeterministic(t) --> false
  end

  context("Purpose: DeterministicDrinkingWater and DeterministicEnvironmentalFlow") do
    dw = DeterministicDrinkingWater(a * 1.m^3/s)
    ef = DeterministicEnvironmentalFlow(b * 1.m^3/s)

    @fact getId(dw) --> not(getId(ef))
    @fact getNeed(dw) --> a * 1e-6
    @fact getNeed(ef) --> b * 1e-6
    @fact isHydropower(dw) --> false
    @fact isHydropower(ef) --> false
    @fact isWaterWithdrawal(dw) --> true
    @fact isWaterWithdrawal(ef) --> true
    @fact isDeterministic(dw) --> true
    @fact isDeterministic(ef) --> true

    @fact getNeed(dw, Day(1)) --> roughly(86400 * a * 1e-6)
    @fact getNeed(ef, Day(1)) --> roughly(86400 * b * 1e-6)
    @fact getNeed(dw, Week(1)) --> roughly(86400 * a * 1e-6 * 7)
    @fact getNeed(ef, Week(1)) --> roughly(86400 * b * 1e-6 * 7)

    @fact need(dw) --> getNeed(dw)
    @fact need(ef) --> getNeed(ef)
  end

  context("Purpose: list constructor (DeterministicPurposes)") do
    dw = DeterministicDrinkingWater(a_discharge)
    ef = DeterministicEnvironmentalFlow(b_discharge)
    list = DeterministicPurposes(drinkingWater=a_discharge, environmentalFlow=b_discharge)
    @fact in(dw, list) --> true
    @fact in(ef, list) --> true
    @fact length(list) --> 2

    list = DeterministicPurposes(drinkingWater=a_discharge)
    @fact in(dw, list) --> true
    @fact in(ef, list) --> false
    @fact length(list) --> 1

    list = DeterministicPurposes(environmentalFlow=b_discharge)
    @fact in(dw, list) --> false
    @fact in(ef, list) --> true
    @fact length(list) --> 1
  end

  context("Purpose: DeterministicHydropower") do
    hp = HydropowerUnit(max_discharge=a_discharge, max_power=b_power, efficiency=1./c)
    p = DeterministicHydropower(hp)

    @fact getId(p) --> not(getId(hp))
    @fact getHydropowerUnit(p) --> hp
    @fact getNeed(p) --> 0.
    @fact getNeed(p, Day(1)) --> 0.
    @fact getNeed(p, Week(1)) --> 0.
    @fact isHydropower(p) --> true

    @fact need(p) --> getNeed(p)
    @fact hydropower(p) --> getHydropowerUnit(p)
  end

  function outputTestShortVersions(g::DamOutput)
    @fact isConditional(g) --> hasCondition(g)
    @fact id(g) --> getId(g)
    @fact condition(g) --> getCondition(g)
  end

  function outputTestShortVersions(g::ConveyanceDamOutput)
    @fact isConditional(g) --> hasCondition(g)
    @fact id(g) --> getId(g)
    @fact condition(g) --> getCondition(g)
    @fact conveyance(g) --> getConveyanceFactor(g)
  end

  function outputTestShortVersions(g::ConstantDamOutput)
    @fact isConditional(g) --> hasCondition(g)
    @fact id(g) --> getId(g)
    @fact condition(g) --> getCondition(g)
    @fact discharge(g) --> getDischarge(g)
  end

  function outputTestShortVersions(g::HydropowerDamOutput)
    @fact isConditional(g) --> hasCondition(g)
    @fact id(g) --> getId(g)
    @fact condition(g) --> getCondition(g)
    @fact hydropower(g) --> getHydropowerUnit(g)
  end

  context("DamOutput: ConstantDamOutput") do
    dout1 = ConstantDamOutput(a_discharge)
    @fact hasCondition(dout1) --> false
    @fact getDischarge(dout1) --> _from_unitful(a_discharge)
    @fact isConstant(dout1) --> true
    @fact isHydropower(dout1) --> false
    @fact hasConveyanceFactor(dout1) --> false
    outputTestShortVersions(dout1)

    cond = MinimumReservoirLevelCondition(a_height)
    dout2 = ConstantDamOutput(a_discharge, cond)
    @fact hasCondition(dout2) --> true
    @fact condition(dout2) --> exactly(cond)
    @fact getDischarge(dout2) --> _from_unitful(a_discharge)
    @fact isConstant(dout2) --> true
    @fact isHydropower(dout2) --> false
    @fact hasConveyanceFactor(dout2) --> false
    outputTestShortVersions(dout2)

    cond = MinimumReservoirLevelCondition(a_height)
    dout3 = ConstantDamOutput(an_id, a_discharge, cond)
    @fact getId(dout3) --> an_id
    @fact hasCondition(dout3) --> true
    @fact condition(dout3) --> exactly(cond)
    @fact getDischarge(dout3) --> _from_unitful(a_discharge)
    @fact isConstant(dout3) --> true
    @fact isHydropower(dout3) --> false
    @fact hasConveyanceFactor(dout3) --> false
    outputTestShortVersions(dout3)
  end

  context("DamOutput: HydropowerDamOutput") do
    hp = HydropowerUnit(max_discharge=a_discharge, max_power=b_power, efficiency=1./c)
    dout1 = HydropowerDamOutput(hp)
    @fact hasCondition(dout1) --> false
    @fact getDischarge(dout1) --> getMaximumDischarge(hp)
    @fact isConstant(dout1) --> false
    @fact isHydropower(dout1) --> true
    @fact hasConveyanceFactor(dout1) --> false
    @fact getId(dout1) --> not(getId(hp))
    outputTestShortVersions(dout1)

    cond = MinimumReservoirLevelCondition(a_height)
    dout2 = HydropowerDamOutput(hp, cond)
    @fact hasCondition(dout2) --> true
    @fact condition(dout2) --> exactly(cond)
    @fact getDischarge(dout2) --> _from_unitful(a_discharge)
    @fact isConstant(dout2) --> false
    @fact isHydropower(dout2) --> true
    @fact hasConveyanceFactor(dout2) --> false
    @fact getId(dout2) --> not(getId(hp))
    outputTestShortVersions(dout2)
  end

  context("DamOutput: ConveyanceSpillway") do
    dout1 = ConveyanceSpillway(a)
    @fact hasCondition(dout1) --> false
    @fact getConveyanceFactor(dout1) --> a
    @fact isConstant(dout1) --> false
    @fact isHydropower(dout1) --> false
    @fact hasConveyanceFactor(dout1) --> true
    outputTestShortVersions(dout1)

    cond = MinimumReservoirLevelCondition(a_height)
    dout2 = ConveyanceSpillway(a, cond)
    @fact hasCondition(dout2) --> true
    @fact condition(dout2) --> exactly(cond)
    @fact getConveyanceFactor(dout2) --> a
    @fact isConstant(dout2) --> false
    @fact isHydropower(dout2) --> false
    @fact hasConveyanceFactor(dout2) --> true
    outputTestShortVersions(dout2)

    cond = MinimumReservoirLevelCondition(a_height)
    dout3 = ConveyanceSpillway(an_id, a, cond)
    @fact getId(dout3) --> an_id
    @fact hasCondition(dout3) --> true
    @fact condition(dout3) --> exactly(cond)
    @fact getConveyanceFactor(dout3) --> a
    @fact isConstant(dout3) --> false
    @fact isHydropower(dout3) --> false
    @fact hasConveyanceFactor(dout3) --> true
    outputTestShortVersions(dout3)
  end

  context("DamOutput: general implementation of methods for user-defined releases") do
    # Use the existing method, based on the id field.
    t1 = Test1DamOutput(an_id)
    @fact getId(t1) --> an_id
    @fact id(t1) --> an_id

    # Use the user-defined getId.
    t1 = Test2DamOutput(an_id)
    @fact getId(t1) --> an_id
    @fact id(t1) --> an_id
  end
end
