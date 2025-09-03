
function class_blk(dh, blk, ct, advsearch, class_rule)
  rangedist = ustrip(advsearch.ball.radii[1])
  p = centroid(blk)
  n = search(p, advsearch)
  if length(n) == advsearch.k
    dists = [evaluate(Euclidean(), ustrip.(to(p)), ustrip.(to(centroid(dh, ind)))) for ind in n]
    class_rule(dists, rangedist) ? ct : 4
  else
    4
  end
end
