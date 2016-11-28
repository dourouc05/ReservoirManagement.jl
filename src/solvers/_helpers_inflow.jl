function hlp_JuMP_inflow_discharge(model::JuMP.Model, r::Reservoir, s::Int, t::Int;
                                   imposeDiverted::Union{Bool, JuMP.JuMPArray}=true)
  if imposeDiverted === false # Comparison on value and type.
    error("When not imposing the value of the diverted rivers, must provide a set of variables.")
  end

  return sum([hlp_JuMP_inflow_discharge(model, r, riv, s, t, imposeDiverted=imposeDiverted) for riv in eachriver(r)])
end

hlp_JuMP_inflow_discharge(model::JuMP.Model, r::Reservoir, nr::NaturalRiver, s::Int, t::Int; kwargs...) =
  getScenario(nr, s, t)

function hlp_JuMP_inflow_discharge(model::JuMP.Model, r::Reservoir, dr::DivertedRiver, s::Int, t::Int;
                                   imposeDiverted::Union{Bool, JuMP.JuMPArray, JuMP.Variable}=true)
  if imposeDiverted === false # Comparison on value and type.
    error("When not imposing the value of the diverted rivers, must provide a set of variables in the argument imposeDiverted.")
  end

  if imposeDiverted === true # Comparison on value and type.
    return maxAllowableFlow(dr, s, t)
  else
    if isa(imposeDiverted, JuMP.JuMPArray)
      imposeDiverted = imposeDiverted[s, t, getName(dr)]
    end

    @constraint(model, imposeDiverted <= maxAllowableFlow(dr, s, t))
    return imposeDiverted
  end
end
