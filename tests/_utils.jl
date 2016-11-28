gen_scenario(length::Int) = Float64[2 + sin(2π * t / 52) + rand() for t in 1:length]
gen_dates(length::Int, year::Int) = DateTime[DateTime(year, October, 1) + Week(i) for i in 0:length - 1]
gen_ts(length::Int, year::Int) = TimeArray(gen_dates(length, year), gen_scenario(length), ["Discharge (m³/s)"], Week(1))
gen_tses(length::Int, year::Int, n::Int) = [gen_ts(length, year + i - 1) for i in 1:n]
gen_river(name::AbstractString, length::Int, year::Int, n::Int) =
  NaturalRiver(name=name, scenarios=gen_tses(length, year, n))
gen_river(name::AbstractString, length::Int, year::Int, n::Int, envFlow::Discharge, maxFlow::Discharge) =
  DivertedRiver(name=name, scenarios=gen_tses(length, year, n), environmental_flow=envFlow, maximum_flow=maxFlow)

function compare_rivers(r1::NaturalRiver, r2::NaturalRiver)
  @fact getName(r1) --> getName(r2)
  @fact countScenarios(r1) --> countScenarios(r2)
  @fact countTimeSteps(r1) --> countTimeSteps(r2)
  @fact getPeriod(r1) --> getPeriod(r2)

  for s in 1:countScenarios(r1)
    @fact timestamp(getScenarioTS(r1, s)) --> timestamp(getScenarioTS(r2, s))
    @fact getScenario(r1, s) --> getScenario(r2, s)
    @fact colnames(getScenarioTS(r1, s)) --> colnames(getScenarioTS(r2, s))
    @fact meta(getScenarioTS(r1, s)) --> meta(getScenarioTS(r2, s))
  end
end
