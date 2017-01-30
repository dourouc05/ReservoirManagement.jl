# Scenario generation

## High-level interface

```@docs
scenarioGeneration(::Reservoir, ::Union{Symbol, Function}, number::Int)
scenarioGeneration(::River, ::Symbol, number::Int)
```

The following scenario generation techniques are available by default:

  * `:Duplicate`: each year is duplicated `number` times
  * `:Concatenate`: `number` successive years are concatenated (no scenario is ever used more than once)
  * `:Merge`: `number` successive years are concatenated (every possible combination is returned)
  * `:Mix`: each output scenario is a combination of `number` scenarios; all possible combinations are returned (some scenarios will contain only copies of one scenario, others consecutive scenarios, etc.)

For example, if there are three scenarios (`s1`, `s2`, and `s3`) and if `number` is `2`, the outputs will be:

  * `:Duplicate`: `s1 s1`, `s2 s2`, and `s3 s3`
  * `:Concatenate`: `s1 s2`
  * `:Merge`: `s1 s2` and `s2 s3`
  * `:Mix`: `s1 s1`, `s1 s2`, `s1 s3`, `s2 s1`, `s2 s2`, `s2 s3`, `s3 s1`, `s3 s2`, and `s3 s3`

## Low-level interface

```@docs
scenarioDuplication(::Reservoir, ::Int)
scenarioDuplication(::River, ::Int)
scenarioConcatenation(::Reservoir, ::Int)
scenarioConcatenation(::River, ::Int)
scenarioMerging(::Reservoir, ::Int)
scenarioMerging(::River, ::Int)
scenarioMixing(::Reservoir, ::Int)
scenarioMixing(::River, ::Int)
```

## Scenario generation technique implementation

Each of the four implemented technique is implemented as a function that takes as input an array of time series (on series being a scenario) and a number, and then returns an array of time series according to its generation algorithm. Each new technique must be registered with the function `addScenarioGenerationAlgorithm!`, that associates a symbol (like `:Duplicate`, `:Mix`, etc.), the previously implemented function, and a string (a name for the algorithm).

```@docs
scenarioDuplication(::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, ::Int)
scenarioConcatenation(::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, ::Int)
scenarioMerging(::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, ::Int)
scenarioMixing(::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, ::Int)
addScenarioGenerationAlgorithm!(::Symbol, ::Function, ::AbstractString)
```

## Utility functions

```@docs
scenarioShift(::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1})
scenarioShift(::River)
scenarioShift(::Reservoir)
```
