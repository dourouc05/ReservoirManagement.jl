"""
Describes the reservoir bathymetry, i.e. the relationship between the reservoir level and the total storage.
The function is only valid within the bounds of the reservoir (i.e. between its minimum and maximum volumes).
The parameter is the volume and its output is the corresponding reservoir level. The function must be one-way.

`isempty(bathymetry)` allows checking whether the bathymetry is actually available, with `NullReservoirBathymetry` being
a nonexistent bathymetry (i.e. no information available).

`computeLevel(bathymetry, v)` gets the reservoir level from the given volume `v` (10^6 m^3).
`computeVolume(bathymetry, h)` get the volume from the given reservoir level `h` (m).
"""
abstract ReservoirBathymetry

immutable NullReservoirBathymetry <: ReservoirBathymetry end
Base.isempty(g::ReservoirBathymetry) = false
Base.isempty(g::NullReservoirBathymetry) = true

computeVolume(g::ReservoirBathymetry, h::Float64) =
  error("Operation not implemented: volume = f(reservoir level), for bathymetry " * string(typeof(g)))
volume(g::ReservoirBathymetry, h::Float64) = computeVolume(g, h)
computeLevel(g::ReservoirBathymetry, v::Float64) =
  error("Operation not implemented: volume = f(reservoir level), for bathymetry " * string(typeof(g)))
level(g::ReservoirBathymetry, h::Float64) = computeLevel(g, h)

computeVolume(g::ReservoirBathymetry, h::Height) = computeLevel(g, _from_unitful(h)) * 1.m^3
computeLevel(g::ReservoirBathymetry, v::Volume) = computeLevel(g, _from_unitful(v)) * 1.m



"""
This special kind of relationship is completely transparent, and can directly be used by an optimisation solver
(if the mathematical structure is allowable).
"""
abstract WhiteBoxReservoirBathymetry <: ReservoirBathymetry

iswhitebox(g::ReservoirBathymetry) = false
iswhitebox(g::WhiteBoxReservoirBathymetry) = true
isblackbox(g::ReservoirBathymetry) = ! iswhitebox(g)



"""
Among the white-box geometries are polynomial-based ones, which can be represented as follows:
    level(volume) = \sum_{i=0}^N a_{i+1} volume^i

These geometries are thus completely defined by a series of coefficients, which *should* be stored
in the `coefficients` field with increasing order (publicly accessible with the `coefficient` function).
In this case, the functions `coefficients` and `level` are automatically provided (not `volume`, as it requires
solving a polynomial equation).
"""
abstract PolynomialReservoirBathymetry <: WhiteBoxReservoirBathymetry

ispolynomial(g::ReservoirBathymetry) = false
ispolynomial(g::PolynomialReservoirBathymetry) = true

isPolynomial(g::ReservoirBathymetry) = ispolynomial(g)

"""
Retrieves the coefficient of the polynomial in `g` at the given order, i.e. power of the volume.
For example, the constant part of the polynomial has order zero, and the linear one has order one.
"""
coefficient(g::PolynomialReservoirBathymetry, order::Int) = g.coefficients[order + 1]

"""
Evaluates the given polynomial reservoir bathymetry at the given volume(s).
"""
function computeLevel(g::PolynomialReservoirBathymetry, v::Float64)
  if length(g.coefficients) == 0
    return zero(volume)
  end

  # Inefficient implementation of polynomial evaluation. Works for single numbers and vectors.
  level = g.coefficients[1]
  for i = 2:length(g.coefficients)
      level += g.coefficients[i] * v .^ (i - 1)
  end
  return level
end



immutable LinearReservoirBathymetry <: PolynomialReservoirBathymetry
  coefficients::Tuple{Float64, Float64}
  LinearReservoirBathymetry(c::Tuple{Float64, Float64}) = new(c)
  LinearReservoirBathymetry(a::Float64, b::Float64) = new((a, b))
end

getIntercept(g::LinearReservoirBathymetry) = g.coefficients[1]
getSlope(g::LinearReservoirBathymetry) = g.coefficients[2]
islinear(g::ReservoirBathymetry) = false
islinear(g::LinearReservoirBathymetry) = true

computeVolume(g::LinearReservoirBathymetry, h::Float64) = (h - getIntercept(g)) / getSlope(g)

isLinear(g::ReservoirBathymetry) = islinear(g)
intercept(g::LinearReservoirBathymetry) = getIntercept(g)
slope(g::LinearReservoirBathymetry) = getSlope(g)




"""
Quadratic reservoir bathymetry, whose form is:
    level = a + b volume + c volume^2
Both convex and concave polynomials make sense, depending on the position of the inflexion point with respect to the
part of the domain that is used:
        level                     level
          ^                         ^
          |    /                    |  ___
          |   /                     | /
          |_ /                      |/
          |                         |
          +------> vol              +------> vol
          Convex polynomial         Concave polynomial
It is the user's responsibility to ensure that this polynomial makes sense, i.e. is strictly increasing between the
reservoir capacities
"""
immutable QuadraticReservoirBathymetry <: PolynomialReservoirBathymetry
  coefficients::Tuple{Float64, Float64, Float64}
  QuadraticReservoirBathymetry(c::Tuple{Float64, Float64, Float64}) = new(c)
  QuadraticReservoirBathymetry(a::Float64, b::Float64, c::Float64) = new((a, b, c))
end

isquadratic(g::ReservoirBathymetry) = false
isquadratic(g::QuadraticReservoirBathymetry) = true

function computeVolume(g::QuadraticReservoirBathymetry, h::Float64)
  # Solve a quadratic equation:
  #    a v^2 + b v + c == h
  # But the well-known formulae only work for == 0, hence the c.
  a = coefficient(g, 2)
  b = coefficient(g, 1)
  c = coefficient(g, 0) - h
  Δ = b^2 - 4 * a * c

  if Δ < 0.
    error("The bathymetry has no real root for h = " * string(h) * ".")
  end

  v1 = (- b + sqrt(Δ)) / (2 * a)
  v2 = (- b - sqrt(Δ)) / (2 * a)

  if v1 < 0. && v2 < 0.
    error("The bathymetry has no acceptable root for h = " * string(h) * ".")
  elseif v1 != v2 && v1 > 0. && v2 > 0.
    error("The bathymetry has two acceptable roots for h = " * string(h) * ".")
  elseif (v1 == 0. && v2 <= 0.) || (v2 == 0. && v1 <= 0.)
    error("The bathymetry has one root for h = " * string(h) * ", but it is zero.")
  elseif v1 == v2
    return v1
  elseif v1 > 0. && v2 <= 0.
    return v1
  elseif v2 > 0. && v1 <= 0.
    return v2
  end
end

isQuadratic(g::ReservoirBathymetry) = isquadratic(g)



immutable PiecewiseLinearReservoirBathymetry <: WhiteBoxReservoirBathymetry
  control_points::Array{Float64, 1} # Strictly increasing.
  values::Array{Float64, 1} # Strictly increasing (reservoir level increases with volume of water).

  function PiecewiseLinearReservoirBathymetry(control_points::Array{Float64, 1}, values::Array{Float64, 1})
    if length(control_points) < 2
      error("Must provide at least a segment, i.e. two control points.")
    end

    if length(control_points) - 1 > length(values)
      error("Inconsistent set of points: too many control points for the values (exactly one more than the values is needed).")
    elseif length(control_points) - 1 < length(values)
      error("Inconsistent set of points: too few control points for the values (exactly one more than the values is needed).")
    end

    for i in 2:length(control_points)
      if control_points[i] <= control_points[i - 1]
        error("Control points must be strictly increasing: problem found near index " * string(i) * ".")
      end
    end

    for i in 2:length(values)
      if values[i] <= values[i - 1]
        error("Values must be strictly increasing: problem found near index " * string(i) * ".")
      end
    end

    return new(control_points, values)
  end
end

ispiecewiselinear(g::ReservoirBathymetry) = false
ispiecewiselinear(g::PiecewiseLinearReservoirBathymetry) = true

isPiecewiseLinear(g::ReservoirBathymetry) = ispiecewiselinear(g)
isPWL(g::ReservoirBathymetry) = ispiecewiselinear(g)



"""
Among the white-box geometries is the square root, which can be represented as follows:
    level(volume) = a_0 + \sqrt{ a_1 \times volume }
As a consequence, the volume is given by:
    volume(level) = \frac{(level - a_0)^2}{a_1}
"""
immutable SquareRootReservoirBathymetry <: WhiteBoxReservoirBathymetry
  a0::Float64
  a1::Float64
  SquareRootReservoirBathymetry(c::Tuple{Float64, Float64}) = new(c[1], c[2])
  SquareRootReservoirBathymetry(a::Float64, b::Float64) = new(a, b)
end

issquareroot(g::ReservoirBathymetry) = false
issquareroot(g::SquareRootReservoirBathymetry) = true

isSquareRoot(g::ReservoirBathymetry) = issquareroot(g)

computeLevel(g::SquareRootReservoirBathymetry, v::Float64) = g.a0 + sqrt(g.a1 * v)
computeVolume(g::SquareRootReservoirBathymetry, h::Float64) = (h - g.a0)^2 / g.a1
