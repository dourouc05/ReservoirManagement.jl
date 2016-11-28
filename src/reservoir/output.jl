"""
Describes the maximum discharge for a dam output, which can be constant (`Constant*DamOutput`) or depend on the
hydraulic head with a conveyance factor (`Conveyance*DamOutput`); the latter make the assumption that the reservoir level
is the hydraulic head. This discharge can be conditional or not, which is indicated by the `condition` attribute
(whose type must be `DamOutputCondition`).

These releases are considered as degrees of freedom of the operator, including for flood management: the operator
can decide on the exact discharge through each dam output (if the conditions are met). These releases are not
supposed to fulfil any dam purpose (like drinking water), which are either constant or valued (see `Purpose` for these).

Each output must have a unique identifier, whose type is `Symbol`. It is always available with the `getId()` function.
For implementers, a default implementation uses the field `id`. Likewise, the condition is available through
`getCondition()`, and a default implementation retrieves the value from the field `condition`.
"""
abstract DamOutput

getId(dout::DamOutput) = dout.id
getCondition(dout::DamOutput) = dout.condition
hasCondition(dout::DamOutput) = ! isempty(getCondition(dout))

id(dout::DamOutput) = getId(dout)
condition(dout::DamOutput) = getCondition(dout)
isConditional(dout::DamOutput) = hasCondition(dout)

"""
Describes the maximum discharge for a dam output using a conveyance factor, be it an open-channel spillway
(`ConveyanceSpillway`) or a closed conduit (`ConveyancePipe`).

A default implementation of `getConveyanceFactor()` takes its value from the `conveyance` field.
"""
abstract ConveyanceDamOutput <: DamOutput

hasConveyanceFactor(dout::DamOutput) = false
hasConveyanceFactor(dout::ConveyanceDamOutput) = true
getConveyanceFactor(dout::ConveyanceDamOutput) = dout.conveyance
conveyance(dout::DamOutput) = getConveyanceFactor(dout)

"""
Condition at which a dam output is enabled. For example, a sufficient reservoir level is needed for the spillways to be
reachable.
"""
abstract DamOutputCondition

"""
Describes a dam output condition: the dam release is always possible.
"""
immutable NoCondition <: DamOutputCondition; end

isempty(c::DamOutputCondition) = false
isempty(c::NoCondition) = true



"""
Describes a constant dam output, independent of all parameters. Its discharge is given in 10^6 m^3/s.
"""
immutable ConstantDamOutput <: DamOutput
  id::Symbol
  discharge::Float64
  condition::DamOutputCondition

  ConstantDamOutput(discharge::Discharge, condition::DamOutputCondition=NoCondition()) =
    new(:constant, _from_unitful(discharge), condition)
  ConstantDamOutput(id::Symbol, discharge::Discharge, condition::DamOutputCondition=NoCondition()) =
    new(id, _from_unitful(discharge), condition)
end

getDischarge(dout::DamOutput) = error("Discharge has no meaning for output " * string(getId(dout)) * ".")
getDischarge(dout::ConstantDamOutput) = dout.discharge
getDischarge(dout::DamOutput, p::Period) = getDischarge(dout) * Dates.days(p) * 86400
isConstant(dout::DamOutput) = false
isConstant(dout::ConstantDamOutput) = true

discharge(dout::DamOutput) = getDischarge(dout)
discharge(dout::DamOutput, p::Period) = getDischarge(dout, p)

"""
Describes the penstock driving a hydropower unit.

This object is automatically created when creating a `Reservoir` data structure, and thus does not need to be used
outside internal code.
"""
immutable HydropowerDamOutput <: DamOutput
  hp::HydropowerUnit
  condition::DamOutputCondition

  HydropowerDamOutput(hp::HydropowerUnit, cond::DamOutputCondition=NoCondition()) = new(hp, cond)
end

getId(hp::HydropowerDamOutput) = getDamOutputId(hp.hp)
getDischarge(hp::HydropowerDamOutput) = getMaximumDischarge(hp.hp)

isHydropower(dout::DamOutput) = false
isHydropower(hp::HydropowerDamOutput) = true
getHydropowerUnit(hp::HydropowerDamOutput) = hp.hp

hydropower(hp::HydropowerDamOutput) = getHydropowerUnit(hp)

"""
Describes the maximum flow of an open-channel spillway using a conveyance factor. This discharge is limited by the
available hydraulic head.

This kind of dam output obeys the following law:
    Q = conveyance (head)^(3/2)

TODO: How to deal with units? Q is 10^6 m^3/period!
"""
immutable ConveyanceSpillway <: ConveyanceDamOutput
  id::Symbol
  conveyance::Float64
  condition::DamOutputCondition

  ConveyanceSpillway(conveyance::Float64, condition::DamOutputCondition=NoCondition()) =
    new(:spillway, conveyance, condition)
  ConveyanceSpillway(id::Symbol, conveyance::Float64, condition::DamOutputCondition=NoCondition()) =
    new(id, conveyance, condition)
end

"""
Describes the maximum flow of a closed conduit using a conveyance factor. This discharge is limited by the available
hydraulic head.

This kind of dam output obeys the following law:
    Q = conveyance (head)^(1/2)

TODO: Is this law the right one?
TODO: How to deal with units? Q is 10^6 m^3/period, while the coefficient will be given for m^3/s!
"""
immutable ConveyancePipe <: ConveyanceDamOutput
  id::Symbol
  conveyance::Float64
  condition::DamOutputCondition

  ConveyancePipe(conveyance::Float64, condition::DamOutputCondition=NoCondition()) =
    new(:penstock, conveyance, condition)
  ConveyancePipe(id::Symbol, conveyance::Float64, condition::DamOutputCondition=NoCondition()) =
    new(id, conveyance, condition)
end



"""
Describes a dam output condition: the release is possible only if this condition is met, i.e. if the reservoir level is
sufficiently high. This type represents physical constraints on the actual reservoir level (not an approximation for
hydraulic head).
"""
immutable MinimumReservoirLevelCondition <: DamOutputCondition
  min::Float64
  MinimumReservoirLevelCondition(min::Height) = new(_from_unitful(min))
end

getMinimumLevel(c::MinimumReservoirLevelCondition) = c.min

getMinimumLevel(c::DamOutputCondition) = error("Minimum reservoir level has no meaning for the condition.")
minLevel(c::DamOutputCondition) = getMinimumLevel(c)

# Equivalent functions for volumes: see output_fwd.jl (they need the declarations of other stuff). 



"""
Builds an array of dam outputs for the most common cases, i.e. one set of bottom outlets and one spillway.
"""
function ConstantDamOutputs(;bottomOutlets::Discharge=0.m^3/s,
                             spillway::Discharge=0.m^3/s, spillway_minLevel::Height=0.m)
  list = DamOutput[] # All data structures expect a generic array of dam outputs, be they constant or not.
  if bottomOutlets != 0.m^3/s
    push!(list, ConstantDamOutput(:bottomOutlets, bottomOutlets))
  end
  if spillway != 0.m^3/s
    if spillway_minLevel != 0.m
      push!(list, ConstantDamOutput(:spillway, spillway, MinimumReservoirLevelCondition(spillway_minLevel)))
    else
      push!(list, ConstantDamOutput(:spillway, spillway))
    end
  end
  return list
end
