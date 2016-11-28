"""
Describes a dam purpose, such as providing drinking water (i.e. water withdrawn from the reservoir), or ensuring inland
navigation (i.e. constraints on the reservoir level behind the dam).

Each purpose must have a unique identifier, whose type is `Symbol`. It is always available with the `getId()` function.
For implementers, a default implementation uses the field `id`.
"""
abstract Purpose

getId(p::Purpose) = p.id
isDeterministic(p::Purpose) = false

id(p::Purpose) = getId(p)

"""
Described a dam purpose that is fulfilled by taking water from the reservoir, such as providing drinking water.

The given discharges are expressed as 10^6 m^3/s.
"""
abstract WaterWithdrawalPurpose <: Purpose

isWaterWithdrawal(p::Purpose) = false
isWaterWithdrawal(p::WaterWithdrawalPurpose) = true



"""
Describes a deterministic dam purpose, without any uncertainty, which only withdraws water.

Each such purpose corresponds to a given needed discharge. It is always available with the `getNeed()` function.
For implementers, a default implementation uses the field `need`. The given discharges are expressed as 10^6 m^3/s.

The rationale behind this layer of abstraction (instead of making it concrete, with the various purposes having
a different id) is that some purposes may have multiple instances for the same dam. For example, a conveyance network
may extract water from multiple dams for a given quantity of water to reach; and the same dam may participate in multiple
such summations.
"""
abstract DeterministicWaterWithdrawalPurpose <: WaterWithdrawalPurpose

getNeed(d::DeterministicWaterWithdrawalPurpose) = d.need
getNeed(d::DeterministicWaterWithdrawalPurpose, p::Period) = getNeed(d) * Dates.days(p) * 86400.
isDeterministic(d::DeterministicWaterWithdrawalPurpose) = true

need(d::DeterministicWaterWithdrawalPurpose) = getNeed(d)
need(d::DeterministicWaterWithdrawalPurpose, p::Period) = getNeed(d, p)



"""
Represents a deterministic drinking water consumption from the reservoir. It can only be fulfilled by one reservoir.

By default, the id is `:drinkingWater`.
"""
immutable DeterministicDrinkingWater <: DeterministicWaterWithdrawalPurpose
  id::Symbol
  need::Float64 # 10^6 m^3/s
  DeterministicDrinkingWater(v::Discharge) = new(:drinkingWater, _from_unitful(v))
  DeterministicDrinkingWater(id::Symbol, v::Discharge) = new(id, _from_unitful(v))
end

"""
Represents a deterministic environmental flow that must be released from the reservoir into its river.

By default, the id is `:environmentalFlow`.
"""
immutable DeterministicEnvironmentalFlow <: DeterministicWaterWithdrawalPurpose
  id::Symbol
  need::Float64 # 10^6 m^3/s
  DeterministicEnvironmentalFlow(v::Discharge) = new(:environmentalFlow, _from_unitful(v))
  DeterministicEnvironmentalFlow(id::Symbol, v::Discharge) = new(id, _from_unitful(v))
end

"""
Represents a hydropower unit, which has no need for water as the exact discharge is a degree of freedom.

The id is computed with that of the hydropower unit as `getPurposeId(hp)`.

This object is automatically created when creating a `Reservoir` data structure, and thus does not need to be used
outside internal code.
"""
immutable DeterministicHydropower <: DeterministicWaterWithdrawalPurpose
  hp::HydropowerUnit
end

getId(hp::DeterministicHydropower) = getPurposeId(hp.hp)
getNeed(hp::DeterministicHydropower) = 0.
getHydropowerUnit(hp::DeterministicHydropower) = hp.hp
isHydropower(p::Purpose) = false
isHydropower(p::DeterministicHydropower) = true

hydropower(hp::DeterministicHydropower) = hp.hp



"""
Builds an array of deterministic purposes for the most common cases, i.e. one `DeterministicDrinkingWater`
and one `DeterministicEnvironmentalFlow` at most.
"""
function DeterministicPurposes(;drinkingWater::Discharge=0.m^3/s, environmentalFlow::Discharge=0.m^3/s)
  list = Purpose[] # All data structures expect a generic array of purposes, be they deterministic or not.
  if drinkingWater != 0.m^3/s
    push!(list, DeterministicDrinkingWater(drinkingWater))
  end
  if environmentalFlow != 0.m^3/s
    push!(list, DeterministicEnvironmentalFlow(environmentalFlow))
  end
  return list
end
