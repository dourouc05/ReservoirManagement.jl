"""
Implements the needed discharge for each purpose. (Purposes without discharge do not need to implement this function.)
The returned discharge may be an expression based on other variables.

For example, imposing a constraint on the reservoir level within the reservoir is not a discharge: the corresponding
function should return `0.`.
"""
hlp_JuMP_purposes_discharge(model::JuMP.Model, r::Reservoir, period::Period) =
  hlp_JuMP_purposes_discharge(model, getPurposes(r), period)

hlp_JuMP_purposes_discharge(model::JuMP.Model, purposes::Array{Purpose, 1}, period::Period) =
  sum([hlp_JuMP_purposes_discharge(model, p, period) for p in purposes])

# Generic purposes have no needed discharge.
hlp_JuMP_purposes_discharge(model::JuMP.Model, purpose::Purpose, period::Period) = 0.

hlp_JuMP_purposes_discharge(model::JuMP.Model, purpose::WaterWithdrawalPurpose, period::Period) =
  error("[hlp_JuMP_purposes_discharge] Helper not defined for " * string(typeof(purpose)) * ".")

hlp_JuMP_purposes_discharge(model::JuMP.Model, purpose::DeterministicWaterWithdrawalPurpose, period::Period) =
  getNeed(purpose, period)

# Hydropower has no imposed discharge.
hlp_JuMP_purposes_discharge(model::JuMP.Model, purpose::HydropowerDamOutput, period::Period) = 0.
