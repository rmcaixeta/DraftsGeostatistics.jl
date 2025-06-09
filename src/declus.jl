## Block weighting advanced
oweight(object, method::BlockWeighting, offsets::Int; kwargs...) = oweight(domain(object), method, offsets; kwargs...)

function oweight(domain::Domain, method::BlockWeighting, offsets::Int; mode="simple")
  orig = boundingbox(domain).min
  offd = -1 .* method.sides ./ (offsets + 1)
  weights = zeros(nelements(domain) + 1)
  dims = embeddim(domain)
  offrange = mode == "prod" ? Iterators.product([0:offsets for d in 1:dims]...) : 0:offsets

  for d in offrange
    opt = orig |> Translate((d .* offd)...)
    odomain = vcat(GeometrySet([opt]), domain)
    p = partition(odomain, BlockPartition(method.sides))
    for s in indices(p)
      sn = filter(i -> i != 1, s)
      weights[sn] .+= 1 / length(sn)
    end
  end

  weights = weights[2:end] ./ length(0:offsets)

  GeoWeights(domain, weights)
end

function cell_declus_tests_i(comps::AbstractGeoTable, evar, spacing, offsets, mode)
  w = oweight(comps, BlockWeighting(spacing...), offsets; mode)
  (; :test => string(spacing), Symbol("mean_$evar") => mean(getproperty(comps, evar), w))
end

function cell_declus_tests(comps::AbstractGeoTable, evar, iterator; offsets=3, mode="simple")
  outvar = Symbol("mean_$evar")
  fun = Transducers.Map(spacing -> cell_declus_tests_i(comps, evar, spacing, offsets, mode))
  res = foldxt(vcat, fun, iterator)
  DataFrame(res) |> Sort(outvar)
end
