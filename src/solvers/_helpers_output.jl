_hlp_outputs_allowed_nonlinear_modes =[:None, :LinearApprox, :PiecewiseLinear, :Cone]
function _hlp_outputs_allowed_nonlinear_mode(mode::Symbol)
  if ! in(mode, _hlp_outputs_allowed_nonlinear_modes)
    error("Given nonlinear mode not recognised: " * string(mode) * ". " *
          "Allowed modes are: " * string(_hlp_outputs_allowed_nonlinear_modes) * ".")
  end
end

"""
Adds the necessary storage for the given `output` variable to implement the asked `outputs` for this `r`eservoir.
This function creates the necessary variables and (conic) constraints, with the given way of handling nonlinearity
for these outputs.

The way of encoding the nonlinearity is configurable by the `nonlinear_mode` variable:

  * `:None`: the nonlinear constraints are forbidden; enabling those outputs will raise an error.
  * `:LinearApprox`: implement those constraints with linear constraints (no binary variables for this).
  * `:PiecewiseLinear`: implement those constraints with a piecewise linear function (with binary variables).
  * `:Cone`: implement those constraints as cones (those that are tolerated by JuMP).

TODO: better way than passing the variables as arguments? Use JuMP ext dictionary to store all variables, replace them just by an index (scenario, time)?
"""
function hlp_JuMP_outputs_constraint(model::JuMP.Model, output::JuMP.Variable, volume::JuMP.Variable,
                                     r::Reservoir, outputs::Array{Symbol, 1}=outputsId(r);
                                     nonlinear_mode::Symbol=:None)
  _hlp_outputs_allowed_nonlinear_mode(nonlinear_mode)
  if length(outputs) >= 1
    @constraint(model, output <= sum([hlp_JuMP_output(model, volume, r, out) for out in outputs]))
  end
  # Else: output <= \infty.
end

"""
Returns the maximum potential release for the given `output` (indexed by id) in this `r`eservoir.
This function creates the necessary variables and (conic) constraints, with the given way of handling nonlinearity
for these outputs.
"""
function hlp_JuMP_output(model::JuMP.Model, volume::JuMP.Variable, r::Reservoir, output::Symbol;
                         nonlinear_mode::Symbol=:None)
  _hlp_outputs_allowed_nonlinear_mode(nonlinear_mode)
  if ! hasOutput(r, output)
    error("The reservoir does not have the output " * string(output) * "!")
  end
  return hlp_JuMP_output(model, volume, r, getOutput(r, output), nonlinear_mode=nonlinear_mode)
end



function hlp_JuMP_output(model::JuMP.Model, volume::JuMP.Variable, r::Reservoir, output::DamOutput;
                         nonlinear_mode::Symbol=:None)
  error("[hlp_JuMP_output] Helper not defined for " * string(typeof(output)) * ".")
end

function hlp_JuMP_output(model::JuMP.Model, volume::JuMP.Variable, r::Reservoir, output::ConstantDamOutput;
                         nonlinear_mode::Symbol=:None)
  # Ignore nonlinear_mode: not relevant for these outputs.
  if isConditional(output)
    return getDischarge(output, getPeriod(r)) * hlp_JuMP_output_condition(model, volume, r, getCondition(output))
  else
    return getDischarge(output, getPeriod(r))
  end
end

function hlp_JuMP_output(model::JuMP.Model, volume::JuMP.Variable, r::Reservoir, output::HydropowerDamOutput;
                         nonlinear_mode::Symbol=:None)
  # Ignore nonlinear_mode: not relevant for these outputs.
  if isConditional(output)
    return getDischarge(output, getPeriod(r)) * hlp_JuMP_output_condition(model, volume, r, getCondition(output))
  else
    return getDischarge(output, getPeriod(r))
  end
end





function hlp_JuMP_output_condition(model::JuMP.Model, volume::JuMP.Variable, r::Reservoir, condition::DamOutputCondition)
  error("[hlp_JuMP_output_condition] Helper not defined for " * string(typeof(DamOutputCondition)) * ".")
end

function hlp_JuMP_output_condition(model::JuMP.Model, volume::JuMP.Variable, r::Reservoir, condition::MinimumReservoirLevelCondition)
  if getMaximumCapacity(r) < getMinimumVolume(r, condition)
    warn("An output may not be activated due to its condition: its minimum reservoir volume " *
         "(" * string(getMinimumVolume(condition)) * ", i.e. a water depth of " * getMinimumLevel(condition) * " m) is higher than the reservoir maximum capacity " *
         "(" * string(getMaximumCapacity(r)) * ") and will thus never be reached. ")
  end

  # Abstract model:
  #     canReleaseWater = 1 if reservoirLevel >= min(cond), 0 otherwise.
  # Implemented as (with stored volume instead of reservoir level):
  #     min(cond) * canReleaseWater <= reservoirStorage
  #     reservoirStorage <= min(cond) * (1 - canReleaseWater) + maxVolume * canReleaseWater
  @variable(model, canReleaseWater, Bin)
  @constraint(model, getMinimumVolume(r, condition) * canReleaseWater <= volume)
  @constraint(model, volume <= getMinimumVolume(r, condition) * (1 - canReleaseWater) + maxCapacity(r) * canReleaseWater)
  return canReleaseWater
end
