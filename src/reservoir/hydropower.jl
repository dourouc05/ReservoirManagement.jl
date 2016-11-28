"""
Describes a hydropower unit.

Hydropower is both a purpose and a possible dam output. The optimisation process could ignore the purpose side of
hydropower, and only consider it as a possible (constant) output; but optimising this purpose needs the dedicated
penstock. As a consequence, the
"""
immutable HydropowerUnit
  id::Symbol
  max_discharge::Float64
  max_power::Float64
  efficiency::Float64

  # Private constructor, without units nor defaut values.
  function HydropowerUnit(id::Symbol, max_discharge::Float64, max_power::Float64, efficiency::Float64)
    if efficiency < 0. || efficiency > 1.
      error("Wrong efficiency: " * string(efficiency) * " not between 0.0 and 1.0")
    end

    new(id, max_discharge, max_power, efficiency)
  end
end

# Public constructor, with units and keyword arguments.
HydropowerUnit(;id::Symbol=:hydropower,
                max_discharge::Union{Discharge, Float64}=error("Missing maximum discharge"),
                max_power::Union{Power, Float64}=error("Missing maximum power not specified"),
                efficiency::Float64=error("Missing efficiency")
          ) = HydropowerUnit(id, _from_unitful(max_discharge), _from_unitful(max_power), efficiency)
HydropowerUnit(max_discharge::Union{Discharge, Float64}, max_power::Union{Power, Float64}, efficiency::Float64) =
  HydropowerUnit(max_discharge=max_discharge, max_power=max_power, efficiency=efficiency)

copy(hp::HydropowerUnit; id::Symbol=hp.id,
     max_discharge::Union{Discharge, Float64}=hp.max_discharge, max_power::Union{Power, Float64}=hp.max_power,
     efficiency::Float64=hp.efficiency) =
  # Always call the internal constructor!
  Reservoir(id, _from_unitful(max_discharge), _from_unitful(max_power), efficiency)




getId(hp::HydropowerUnit) = hp.id
getPurposeId(hp::HydropowerUnit) = symbol(getId(hp), "Purpose")
getDamOutputId(hp::HydropowerUnit) = symbol(getId(hp), "DamOutput")
getMaximumDischarge(hp::HydropowerUnit) = hp.max_discharge
getMaximumPower(hp::HydropowerUnit) = hp.max_power
getEfficiency(hp::HydropowerUnit) = hp.efficiency

id(hp::HydropowerUnit) = getId(hp)
purposeId(hp::HydropowerUnit) = getPurposeId(hp)
damOutputId(hp::HydropowerUnit) = getDamOutputId(hp)
maxDischarge(hp::HydropowerUnit) = getMaximumDischarge(hp)
maxPower(hp::HydropowerUnit) = getMaximumPower(hp)
efficiency(hp::HydropowerUnit) = getEfficiency(hp)
