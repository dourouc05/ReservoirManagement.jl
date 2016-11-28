"""
A generic river that flows down into a reservoir. It might be diverted or not.
This river is named `name`, and is described by a sequence of time series
in `scenarios`. The sampling periods for the same series are indicated by the
`period` variable. Each scenario is one-year-long, and the overall sequence
should be ordered along time (i.e., if the first scenario is taken from year Y,
then the second one must be year Y+1).

Extensions thereof can be implemented in the dictionary `ext` (mapping strings
to objects).
"""
abstract River



### Generic getters.

getName(r::River) = r.name
getScenarios(r::River) = r.scenarios
getScenarioTS(r::River, scenario::Int) = r.scenarios[scenario]
getScenario(r::River, scenario::Int) = values(getScenarioTS(r, scenario))
getScenario(r::River, scenario::Int, time_step::Int) = getScenario(r, scenario)[time_step]
# It does not make sense to add a period argument to this, as it would integrate the value of the scenario
# over the given period (which could overlap with the next time steps). The best way to deal with this is to perform
# an actual integration over the right period when building the object, i.e. in the method makeScenarios().

countScenarios(r::River) = length(getScenarios(r))
countScenariosWithConsistency(r::River) = countScenarios(r)

countTimeSteps(r::River) = length(getScenario(r, 1))
countTimeStepsWithConsistency(r::River) = countTimeSteps(r)

getPeriod(r::River) = meta(r.scenarios[1])
getPeriodWithConsistency(r::River) = getPeriod(r)

getTimePoints(r::River) = timestamp(r.scenarios[1])
getTimePointsWithConsistency(r::River) = getTimePoints(r)

name(r::River) = getName(r)
scenarios(r::River) = getScenarios(r)
scenario(r::River, scenario::Int) = getScenario(r, scenario)
scenario(r::River, scenario::Int, time_step::Int) = getScenario(r, scenario, time_step)
period(r::River) = getPeriod(r)

### Constructor helpers.

function _checkConsistency(r::River)
  if countScenarios(r) < 1
    error("[River " * getName(r) * "] No scenarios!")
  end

  for scenario in 1:countScenarios(r)
    if any(isnan(getScenario(r, scenario)))
      error("[River " * getName(r) * "] NaN found in scenario " * string(scenario) * "!")
    end
    if ! (typeof(meta(getScenarioTS(r, scenario))) <: Base.Dates.Period)
      error("[River " * getName(r) * "] Wrong meta for the input time series: must be a Period; " *
            "got a " * string(typeof(meta(scenarios[1]))) * "!")
    end
  end
end
