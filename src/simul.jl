
function nreal_partitions(nreal, npart)
  flat_v = 1:nreal
  base_size, remainder = divrem(nreal, npart)
  v = Vector{Vector{Int}}()

  start_idx = 1
  for i in 1:npart
    part_size = base_size + (i <= remainder ? 1 : 0)
    end_idx = start_idx + part_size - 1
    reals = flat_v[start_idx:end_idx]
    length(reals) > 0 && push!(v, reals)
    start_idx = end_idx + 1
  end
  v
end

function simsvcat(ensembles...)
  domain = ensembles[1].domain
  fetch_ = ensembles[1].fetch
  reals = mapreduce(e -> e.reals, vcat, ensembles)
  Ensemble(domain, reals, fetch_)
end
