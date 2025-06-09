using GeoStats
using CSV
using DataFrames
using Transducers

nthr = Base.Threads.nthreads()
mapxt = Transducers.Map
geovcat(items...) = vcat([x for x in items if x != nothing]...)

function class_blk(dh, to_estim, ct, advsearch)
  grid = to_estim.geometry
  #knn = KNearestSearch(dh, 1)
  rangedist = ustrip(advsearch.ball.radii[1])
  catgs = mapreduce(vcat, grid) do blk
    p = centroid(blk)
    n = search(p, advsearch)
    outval = if length(n) == 3
      dists = [evaluate(Euclidean(), ustrip.(to(p)), ustrip.(to(centroid(dh, ind)))) for ind in n]
      (sum(dists) < rangedist * 2) ? ct : NaN
    else
      NaN
    end
    outval
  end
  catgs = catgs isa AbstractVector ? catgs : [catgs]
  georef((RCATG=catgs, IJK=to_estim.IJK), grid)
end
