getMinimumVolume(r::Reservoir, c::MinimumReservoirLevelCondition) = computeVolume(bathymetry(r), getMinimumLevel(c))
getMinimumVolume(r::Reservoir, c::DamOutputCondition) = error("Minimum reservoir volume has no meaning for the condition.")
minVolume(r::Reservoir, c::MinimumReservoirLevelCondition) = getMinimumVolume(r, c)
