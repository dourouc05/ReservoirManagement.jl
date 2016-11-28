facts("Scenario generation") do
  context("Internal data structure") do
    @fact _scgen_symbol_from_function(scenarioDuplication) --> :Duplicate
    @fact _scgen_symbol_from_function(scenarioConcatenation) --> :Concatenate
    @fact _scgen_symbol_from_function(scenarioMerging) --> :Merge
    @fact _scgen_symbol_from_function(scenarioMixing) --> :Mix

    @fact _scgen_function_from_symbol(:Duplicate) --> exactly(scenarioDuplication)
    @fact _scgen_function_from_symbol(:Concatenate) --> exactly(scenarioConcatenation)
    @fact _scgen_function_from_symbol(:Merge) --> exactly(scenarioMerging)
    @fact _scgen_function_from_symbol(:Mix) --> exactly(scenarioMixing)

    @fact _scgen_string_from_symbol(:Duplicate) --> "duplicate"
    @fact _scgen_string_from_symbol(:Concatenate) --> "concatenate"
    @fact _scgen_string_from_symbol(:Merge) --> "merge"
    @fact _scgen_string_from_symbol(:Mix) --> "mix"

    @fact _scgen_string_from_function(scenarioDuplication) --> "duplicate"
    @fact _scgen_string_from_function(scenarioConcatenation) --> "concatenate"
    @fact _scgen_string_from_function(scenarioMerging) --> "merge"
    @fact _scgen_string_from_function(scenarioMixing) --> "mix"

    @fact _scgen_string_from_any(:Duplicate) --> "duplicate"
    @fact _scgen_string_from_any(:Concatenate) --> "concatenate"
    @fact _scgen_string_from_any(:Merge) --> "merge"
    @fact _scgen_string_from_any(:Mix) --> "mix"

    @fact _scgen_string_from_any(scenarioDuplication) --> "duplicate"
    @fact _scgen_string_from_any(scenarioConcatenation) --> "concatenate"
    @fact _scgen_string_from_any(scenarioMerging) --> "merge"
    @fact _scgen_string_from_any(scenarioMixing) --> "mix"
  end

  context("Duplication") do
    TS = gen_tses(52, 2000, 3)

    @fact_throws scenarioDuplication([gen_ts(length, year) for i in 1:6], 2) # Same year for each scenario.

    dTS = scenarioDuplication(TS, 2)
    @fact length(dTS) --> length(TS)
    for i in 1:length(dTS)
      @fact values(TS[i][1:52]) --> values(dTS[i][1:52])
      @fact values(TS[i][1:52]) --> values(dTS[i][53:104])
    end

    dTS = scenarioDuplication(TS, 3)
    @fact length(dTS) --> length(TS)
    for i in 1:length(dTS)
      @fact values(TS[i][1:52]) --> values(dTS[i][1:52])
      @fact values(TS[i][1:52]) --> values(dTS[i][53:104])
      @fact values(TS[i][1:52]) --> values(dTS[i][105:156])
    end
  end

  context("Duplication (river APIs)") do
    TS = gen_river("R", 52, 2000, 3)

    dTS = scenarioDuplication(TS, 2)
    @fact countScenarios(dTS) --> countScenarios(TS)
    @fact countTimeSteps(dTS) --> 2 * countTimeSteps(TS)
    for i in 1:countScenarios(dTS)
      @fact getScenario(TS, i)[1:52] --> getScenario(dTS, i)[1:52]
      @fact getScenario(TS, i)[1:52] --> getScenario(dTS, i)[53:104]
    end

    compare_rivers(scenarioGeneration(TS, :Duplicate, 2), dTS)
    compare_rivers(scenarioGeneration(TS, scenarioDuplication, 2), dTS)

    dTS = scenarioDuplication(TS, 3)
    @fact countScenarios(dTS) --> countScenarios(TS)
    @fact countTimeSteps(dTS) --> 3 * countTimeSteps(TS)
    for i in 1:countScenarios(dTS)
      @fact getScenario(TS, i)[1:52] --> getScenario(dTS, i)[1:52]
      @fact getScenario(TS, i)[1:52] --> getScenario(dTS, i)[53:104]
      @fact getScenario(TS, i)[1:52] --> getScenario(dTS, i)[105:156]
    end

    compare_rivers(scenarioGeneration(TS, :Duplicate, 3), dTS)
    compare_rivers(scenarioGeneration(TS, scenarioDuplication, 3), dTS)
  end

  context("Concatenation") do
    TS = gen_tses(52, 2000, 6)

    # The same year cannot be present twice in the same generated scenario. Hence create from scratch a series of
    # time series with the same year, and try to generate scenario from them: this input data is inconsistent.
    @fact_throws scenarioConcatenation([gen_ts(length, year) for i in 1:6], 2)

    @fact_throws scenarioConcatenation(TS, 5) # Not an integer number of input scenarios per output scenario.

    cTS = scenarioConcatenation(TS, 2)
    @fact length(cTS) --> 6 / 2
    @fact values(cTS[1][1:52]) --> values(TS[1][1:52])
    @fact values(cTS[1][53:104]) --> values(TS[2][1:52])
    @fact values(cTS[2][1:52]) --> values(TS[3][1:52])
    @fact values(cTS[2][53:104]) --> values(TS[4][1:52])
    @fact values(cTS[3][1:52]) --> values(TS[5][1:52])
    @fact values(cTS[3][53:104]) --> values(TS[6][1:52])
  end

  context("Concatenation (river APIs)") do
    TS = gen_river("R", 52, 2000, 6)

    dTS = scenarioConcatenation(TS, 2)
    @fact countScenarios(dTS) --> countScenarios(TS) / 2
    @fact countTimeSteps(dTS) --> 2 * countTimeSteps(TS)
    @fact getScenario(dTS, 1)[1:52] --> getScenario(TS, 1)[1:52]
    @fact getScenario(dTS, 1)[53:104] --> getScenario(TS, 2)[1:52]
    @fact getScenario(dTS, 2)[1:52] --> getScenario(TS, 3)[1:52]
    @fact getScenario(dTS, 2)[53:104] --> getScenario(TS, 4)[1:52]
    @fact getScenario(dTS, 3)[1:52] --> getScenario(TS, 5)[1:52]
    @fact getScenario(dTS, 3)[53:104] --> getScenario(TS, 6)[1:52]

    compare_rivers(scenarioGeneration(TS, :Concatenate, 2), dTS)
    compare_rivers(scenarioGeneration(TS, scenarioConcatenation, 2), dTS)
  end

  context("Merging") do
    TS = gen_tses(52, 2000, 4)

    # The same year cannot be present twice in the same generated scenario. Hence create from scratch a series of
    # time series with the same year, and try to generate scenario from them: this input data is inconsistent.
    @fact_throws scenarioMerging([gen_ts(length, year) for i in 1:6], 2)

    cTS = scenarioMerging(TS, 2)
    @fact length(cTS) --> 3
    @fact values(cTS[1][1:52]) --> values(TS[1][1:52])
    @fact values(cTS[1][53:104]) --> values(TS[2][1:52])
    @fact values(cTS[2][1:52]) --> values(TS[2][1:52])
    @fact values(cTS[2][53:104]) --> values(TS[3][1:52])
    @fact values(cTS[3][1:52]) --> values(TS[3][1:52])
    @fact values(cTS[3][53:104]) --> values(TS[4][1:52])
  end

  context("Merging (river APIs)") do
    TS = gen_river("R", 52, 2000, 4)

    dTS = scenarioMerging(TS, 2)
    @fact countScenarios(dTS) --> countScenarios(TS) - 1
    @fact countTimeSteps(dTS) --> 2 * countTimeSteps(TS)
    @fact getScenario(dTS, 1)[1:52] --> getScenario(TS, 1)[1:52]
    @fact getScenario(dTS, 1)[53:104] --> getScenario(TS, 2)[1:52]
    @fact getScenario(dTS, 2)[1:52] --> getScenario(TS, 2)[1:52]
    @fact getScenario(dTS, 2)[53:104] --> getScenario(TS, 3)[1:52]
    @fact getScenario(dTS, 3)[1:52] --> getScenario(TS, 3)[1:52]
    @fact getScenario(dTS, 3)[53:104] --> getScenario(TS, 4)[1:52]

    compare_rivers(scenarioGeneration(TS, :Merge, 2), dTS)
    compare_rivers(scenarioGeneration(TS, scenarioMerging, 2), dTS)
  end

  context("Mixing") do
    TS = gen_tses(52, 2000, 3)

    # The same year cannot be present twice in the same generated scenario. Hence create from scratch a series of
    # time series with the same year, and try to generate scenario from them: this input data is inconsistent.
    @fact_throws scenarioMixing([gen_ts(length, year) for i in 1:6], 2)

    cTS = scenarioMixing(TS, 2)
    @fact length(cTS) --> 9 # All pairs among three objects (can take twice the same object!).

    possible_values = Array{Float64, 1}[values(TS[i][1:52]) for i in 1:3]
    for i in 1:9
      @fact in(values(cTS[i][1:52]), possible_values) --> true
      @fact in(values(cTS[i][53:104]), possible_values) --> true
    end
  end

  context("Mixing (river APIs)") do
    TS = gen_river("R", 52, 2000, 3)

    dTS = scenarioMixing(TS, 2)
    @fact countScenarios(dTS) --> 9
    @fact countTimeSteps(dTS) --> 2 * countTimeSteps(TS)
    @fact getScenario(dTS, 1)[1:52] --> getScenario(TS, 1)[1:52]
    possible_values = Array{Float64, 1}[getScenario(TS, i)[1:52] for i in 1:3]
    for i in 1:9
      @fact in(getScenario(dTS, i)[1:52], possible_values) --> true
      @fact in(getScenario(dTS, i)[53:104], possible_values) --> true
    end

    # The different function calls must give the exact same results, as they are synonymous.
    compare_rivers(scenarioGeneration(TS, :Mix, 2), dTS)
    compare_rivers(scenarioGeneration(TS, scenarioMixing, 2), dTS)
  end
end
