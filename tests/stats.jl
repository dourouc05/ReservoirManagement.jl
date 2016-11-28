facts("Statistical models") do
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

  TSesNR1 = getScenarios(nr1)

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

  # context("Utilities: split an array into random subsets") do
  #   @fact splitIntoSets_sizes(8, [1.]) --> [8]
  #   @fact splitIntoSets_sizes(8, [.8, .2]) --> [6, 2]
  #   @fact_throws splitIntoSets_sizes(8, [1.1]) # Relative sizes do not sum up to 1.
  #
  #   tses = gen_tses(52, 2000, 8)
  #   tses_100 = splitIntoSets(tses, [1.])
  #   @fact length(tses_100) --> 1
  #   @fact length(tses_100[1]) --> 8
  #   tses_80_20 = splitIntoSets(tses, [.8, .2])
  #   @fact length(tses_80_20) --> 2
  #   @fact length(tses_80_20[1]) --> 6
  #   @fact length(tses_80_20[2]) --> 2
  # end
  #
  # context("t-Student: model from weekly data") do
  #   cis = computeCI_TStudent(nr1, .9, ci_period=Week(1))
  #
  #   @fact length(cis) --> 3
  #   for i in 1:3 # mean, low, high
  #     @fact getPeriod(nr1) --> meta(cis[i])
  #     @fact length(cis[i]) --> 52
  #     @fact sum(values(cis[i]) .< 0.) --> 0
  #   end
  #
  #   @fact sum(values(cis[2]) .> values(cis[1])) --> 0 # low < mean
  #   @fact sum(values(cis[2]) .> values(cis[3])) --> 0 # low < high
  #   @fact sum(values(cis[1]) .> values(cis[3])) --> 0 # mean < high
  # end
  #
  # context("Adam2015: model from weekly data") do
  #   cis = computeCI_Adam2015(nr1, .9)
  #
  #   @fact length(cis) --> 3
  #   for i in 1:3 # mean, low, high
  #     @fact getPeriod(nr1) --> meta(cis[i])
  #     @fact length(cis[i]) --> 52
  #     @fact sum(values(cis[i]) .< 0.) --> 0
  #   end
  #
  #   # TODO: not working yet!
  #   # @fact sum(values(cis[2]) .> values(cis[1])) --> 0 # low < mean
  #   # @fact sum(values(cis[2]) .> values(cis[3])) --> 0 # low < high
  #   # @fact sum(values(cis[1]) .> values(cis[3])) --> 0 # mean < high
  # end

  context("Adam2015: scenarios from model (low level interface)") do
    scs = generateScenarios_Adam2015(fit_Adam2015(TSesNR1), Week(1), 2000)
    @fact length(scs) --> 2
    dates = scs[1]
    values = scs[2]
    @fact length(dates) --> length(values)
    @fact length(dates) --> 52
  end

  # context("Adam2015: scenarios from model (high level interface, river)") do
  #   scs = collect(generateScenarios(nr1, 5, Week(1), :Adam2015, Year(1)))
  #   @fact length(scs) --> 5
  #   for s in 1:5
  #     @fact getName(scs[s]) --> getName(nr1)
  #     @fact countScenarios(scs[s]) --> 1
  #     @fact getPeriod(scs[s]) --> Week(1)
  #     @fact countTimeSteps(scs[s]) --> 52
  #     return
  #   end
  # end
end
