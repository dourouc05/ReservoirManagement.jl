"""
A full drainage basin, with a series of rivers and of reservoirs. 

TODO: how to handle topology? 
TODO: unused for now! 
"""
immutable DrainageBasin
  name::AbstractString
  variant::AbstractString

  reservoirs::Array{Reservoir, 1}
  
  function DrainageBasin(name::AbstractString, reservoirs::Array{Reservoir, 1})
    if length(reservoirs) < 1
      error("List of reservoirs does not even contain a reservoir.")
    end
    obj = new(name, "vanilla", reservoirs)
    
    # Check consistency of scenarios. 
    countScenariosWithConsistency(obj)
    countTimeStepsWithConsistency(obj)
    
    return obj
  end
  
  function DrainageBasin(name::AbstractString, reservoir::Reservoir)
    obj = new(name, "vanilla", [reservoir])
    
    # Check consistency of scenarios. 
    countScenariosWithConsistency(obj)
    countTimeStepsWithConsistency(obj)
    
    return obj
  end
end


### Generic getters. 

getName(b::DrainageBasin) = b.name
getVariant(b::DrainageBasin) = b.variant
getReservoirs(b::DrainageBasin) = b.reservoirs

# Counting things. 

"""
Counts the number of scenarios within this basin.
"""
countScenarios(b::DrainageBasin) = countScenarios(getReservoirs(b)[1])

"""
Counts the number of scenarios within this basin and ensures the whole data structure is consistent.
"""
function countScenariosWithConsistency(b::DrainageBasin)
  cnt = -1
  
  for reservoir in getReservoirs(b)
    if cnt < 0
      cnt = countScenariosWithConsistency(reservoir)
    elseif cnt != countScenariosWithConsistency(reservoir)
      error("Inconsistent data in basin " * getName(b) * ": different number of scenarios within the basin (while checking reservoir " * getName(reservoir) * ").")
    end
  end
  
  if cnt < 0
    error("Inconsistent data in basin " * getName(b) * ": no time series found in the rivers.")
  end
  
  return cnt 
end

"""
Counts the number of time steps within this basin.
"""
countTimeSteps(b::DrainageBasin) = countTimeSteps(getReservoirs(b)[1])

"""
Counts the number of time steps within this basin and ensures the whole data structure is consistent.
"""
function countTimeStepsWithConsistency(b::DrainageBasin)
  cnt = -1
  
  for reservoir in getReservoirs(b)
    if cnt < 0
      cnt = countTimeStepsWithConsistency(reservoir)
    elseif cnt != countTimeStepsWithConsistency(reservoir)
      error("Inconsistent data in basin " * getName(b) * ": different number of time steps within the basin (while checking reservoir " * getName(reservoir) * ").")
    end
  end
  
  if cnt < 0
    error("Inconsistent data in basin " * getName(b) * ": no time series found in the rivers.")
  end
  
  return cnt 
end

"""
Retrieves the periodicity of this basin data.
"""
getPeriod(b::DrainageBasin) = getPeriod(getReservoirs(b)[1])

"""
Retrieves the periodicity of this basin data and ensures the whole data structure is consistent.
"""
function getPeriodWithConsistency(b::DrainageBasin)
  period = :Unknown
  
  for reservoir in getReservoirs(b)
    if period == :Unknown
      period = getPeriod(reservoir)
    elseif period != getPeriod(reservoir)
      error("Inconsistent data in basin " * getName(b) * ": different period within the basin (while checking reservoir " * getName(reservoir) * ").")
    end
  end
  
  if period == :Unknown
    error("Inconsistent data in basin " * getName(b) * ": no time series found in the rivers.")
  end
  
  return period::Period
end