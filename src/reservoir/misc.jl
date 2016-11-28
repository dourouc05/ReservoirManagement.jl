countScenarios(TS::Array{TimeSeries.TimeArray{Float64,1,DateTime,Array{Float64,1}}, 1}) = length(TS)
function countTimeSteps(TS::Array{TimeSeries.TimeArray{Float64,1,DateTime,Array{Float64,1}}, 1})
  if countScenarios(TS) == 0
    error("No scenarios are available!")
  end
  return length(TS[1])
end
