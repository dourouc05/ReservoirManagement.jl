"""
Computes the sizes for the actual set splitting in `splitIntoSets`.
"""
function splitIntoSets_sizes(n_scenarios::Int, sizes::Array{Float64, 1})
  @assert abs(sum(sizes) - 1.) < 1e-8 "Relative sizes for each set do not sum up to 1."

  # First approximation of the number of items per set: just round.
  n_each_sets = Array{Int}(round(sizes * n_scenarios))

  # With potential rounding issues, ensure all samples are taken! (No more, no fewer.)
  if sum(n_each_sets) != n_scenarios
    if sum(n_each_sets) < n_scenarios
      # Too *few* scenarios distributed! If some of them are at zero, prefer adding some to them.
      candidates = find(n_each_sets .== 0)
      if length(candidates) > 0
        i = 1
        while sum(n_each_sets) < n_scenarios && i <= length(candidates)
          n_each_sets[candidates[i]] += 1
          i += 1
        end
      end

      # Loop and add elements to each set until convergence.
      i = 1
      while sum(n_each_sets) < n_scenarios
        n_each_sets[i] += 1
        i += 1
      end
    else
      # Too *many* scenarios distributed! Remove elements from the sets that have most of them first.
      candidates = sortperm(n_each_sets)
      i = 1
      while sum(n_each_sets) > n_scenarios
        n_each_sets[candidates[i]] += 1
        i += 1
      end
    end
  end

  @assert sum(n_each_sets) == n_scenarios "Not all scenarios are split in sets!"

  return n_each_sets
end

"""
Splits the various scenarios in `TS` in different sets (e.g., learning set and test set), whose relative sizes are given in `sizes`
(i.e. if `sizes = [.8, .2]`, the output will be two sets, one having 80% of the scenarios, the other the remaining 20%).
The actual separation in sets is random. In other words, this function implements random subsampling.
"""
function splitIntoSets(TS::Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}, sizes::Array{Float64, 1}; forced_perm::Array{Int, 1}=zeros(Int, 0))
  n_scenarios = length(TS)
  n_each_sets = splitIntoSets_sizes(n_scenarios, sizes)

  @assert n_scenarios >= 2 "Not enough scenarios to split in different sets (only one scenario given)."

  # Actually perform the sets separation with a *permutation* of the input scenarios (this is where the randomness comes from).
  if isempty(forced_perm)
    permutation = randperm(n_scenarios)
  else
    permutation = forced_perm
  end
  sc_permuted = TS[permutation]

  cumn = [0; cumsum(n_each_sets)]
  return Array{TimeSeries.TimeArray{Float64, 1, DateTime, Array{Float64, 1}}, 1}[sc_permuted[cumn[i] + 1:cumn[i + 1]] for i in 1:length(n_each_sets)]
end

splitIntoSets(r::NaturalRiver,  sizes::Array{Float64, 1}; kwargs...) =  NaturalRiver[copy(r, scenarios=sc) for sc in splitIntoSets(r.scenarios, sizes; kwargs...)]
splitIntoSets(r::DivertedRiver, sizes::Array{Float64, 1}; kwargs...) = DivertedRiver[copy(r, scenarios=sc) for sc in splitIntoSets(r.scenarios, sizes; kwargs...)]

function splitIntoSets(reservoir::Reservoir, sizes::Array{Float64, 1})
  permutation = randperm(countScenarios(reservoir))
  lnr = [splitIntoSets(r, sizes, forced_perm=permutation) for r in  naturalRivers(reservoir)]
  ldr = [splitIntoSets(r, sizes, forced_perm=permutation) for r in divertedRivers(reservoir)]
  return Reservoir[copy(reservoir, rivers_in=NaturalRiver[lnr[j][i] for j in 1:length(lnr)], rivers_diverted=DivertedRiver[ldr[j][i] for j in 1:length(ldr)]) for i in 1:length(sizes)]
end
