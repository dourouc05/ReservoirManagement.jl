_hlp_bathymetry_allowed_nonlinear_modes =[:None, :PiecewiseLinear, :Cone, :Nonconvex]
function _hlp_bathymetry_allowed_nonlinear_mode(mode::Symbol)
  if ! in(mode, _hlp_outputs_allowed_nonlinear_modes)
    error("Given nonlinear mode not recognised: " * string(mode) * ". " *
          "Allowed modes are: " * string(_hlp_bathymetry_allowed_nonlinear_modes) * ".")
  end
end

"""
Implements the given bathymetry relationship between the reservoir level and the volume in the reservoir.
"""
function hlp_JuMP_bathymetry(model::JuMP.Model, volume::JuMP.Variable, level::JuMP.Variable, r::Reservoir)
  g = bathymetry(r)
  if isempty(g)
    error("The reservoir has no bathymetry!")
  end
  return hlp_JuMP_bathymetry(model, volume, level, bathymetry(r))
end

function hlp_JuMP_bathymetry(model::JuMP.Model, volume::JuMP.Variable, level::JuMP.Variable, g::ReservoirBathymetry)
  error("[hlp_JuMP_bathymetry] Helper not defined for " * string(typeof(g)) * ".")
end

function hlp_JuMP_bathymetry(model::JuMP.Model, volume::JuMP.Variable, level::JuMP.Variable, g::LinearReservoirBathymetry;
                           nonlinear_mode::Symbol=:None)
  # Ignore nonlinear_mode: not relevant for this bathymetry (purely linear).
  @constraint(model, level == intercept(g) + slope(g) * volume)
end

function hlp_JuMP_bathymetry(model::JuMP.Model, volume::JuMP.Variable, level::JuMP.Variable, g::LinearReservoirBathymetry;
                           nonlinear_mode::Symbol=:None)
  _hlp_bathymetry_allowed_nonlinear_mode(nonlinear_mode)

  if mode == :Nonconvex
    @NLconstraint(m, level == coefficient(g, 0) + coefficient(g, 1) * volume + coefficient(g, 2) * storage ^ 2)
  else
    error("[hlp_JuMP_bathymetry] TODO: mode " * string(nonlinear_mode) * " not implemented for LinearReservoirBathymetry.")
  end
end
