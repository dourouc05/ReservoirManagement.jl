# Quick notes for tests. None is actually implemented. This file is not even supposed to work.

using FactCheck

facts("Module: plot") do
  context() do
    # a: no sc, i.e. only one scenario per reservoir.     => with sc, must throw an error.
    # d: no sc, check the second output argument is empty.
    a = getDataFrameForRiver([VesdreReservoir_adam2015_95_avg, VesdreReservoir_adam2015_95_low, VesdreReservoir_adam2015_95_high], 1:ruleLength, "Vesdre", multipleScenarios=:None)
    b = getDataFrameForRiver([sc, VesdreReservoir_adam2015_95_avg, VesdreReservoir_adam2015_95_low, VesdreReservoir_adam2015_95_high], 1:ruleLength, "Vesdre", multipleScenarios=:Merge)
    c = getDataFrameForRiver([sc, VesdreReservoir_adam2015_95_avg, VesdreReservoir_adam2015_95_low, VesdreReservoir_adam2015_95_high], 1:ruleLength, "Vesdre", multipleScenarios=:Split)
    d = getDataFrameForRiver([VesdreReservoir_adam2015_95_avg, VesdreReservoir_adam2015_95_low, VesdreReservoir_adam2015_95_high], 1:ruleLength, "Vesdre", multipleScenarios=:Split)
  end
end

facts("Module: read") do
  context() do
    b = makeScenarios(Vesdre_ts, sampleDuration=Week(1))
    c = makeScenarios(Vesdre_ts, sampleDuration=Day(1), resamplingRate=28) ---> throws()
  end
end

facts("Module: reservoir") do
  context("accessors: statistical") do
    splitIntoSets(Vesdre.scenarios, [1])
    splitIntoSets(Vesdre.scenarios, [.8, .2])
  end

  context("accessors") do
    totalInflow(VesdreReservoir, 1)
    totalInflow(VesdreReservoir, 1:22)
    totalInflow(VesdreReservoir, 1, asTimeSeries=true)
    totalInflow(VesdreReservoir, 1:22, asTimeSeries=true)
  end
end
