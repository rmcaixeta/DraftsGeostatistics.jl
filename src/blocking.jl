
using GeoStats
using Transducers
using DataFrames
using StaticArrays

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

function blocks_iterator(; x=nothing, y=nothing, xy=nothing, z=nothing)
  is2d = isnothing(z)
  issquare = !isnothing(xy)
  isrect = !isnothing(x) && !isnothing(y)
  !issquare && !isrect && error("invalid input")

  ranges = is2d && issquare ? (xy,) : is2d ? (x, y) : issquare ? (xy, z) : (x, y, z)
  iterators = Iterators.product(ranges...)

  mapreduce(vcat, iterators) do i
    #println(i)
    is2d && issquare ? (i[1], i[1]) : is2d ? (i[1], i[2]) : issquare ? (i[1], i[1], i[2]) : (i[1], i[2], i[3])
  end
end

## SUBBLOCK INTERFACE

function get_spacing(obj; coords=[:x, :y, :z])
  spacing = if obj isa CartesianGrid
    obj.spacing
  elseif obj isa AbstractDataFrame
    dx = mapreduce(vcat, coords) do c
      vals = sort(unique(obj[!, c]))
      vals = [vals[i] - vals[i - 1] for i in 2:length(vals)]
      minimum(vals)
    end
    tuple(dx...)
  else
    dom = domain(obj)
    dom = dom isa SubDomain ? dom.domain : dom
    if dom isa CartesianGrid
      dom.spacing
    else
      blk = dom[1]
      verts = to.(blk.vertices)
      dx = mapreduce(vcat, 1:length(verts[1])) do i
        cx = [c[i] for c in verts]
        maximum(cx) - minimum(cx)
      end
      tuple(dx...)
    end
  end
end

function reblock(df::AbstractDataFrame, vars::AbstractVector)
  rvars = [v => mean => v for v in vars]
  combine(groupby(df, :IJK), rvars...)
end

function reblock(gt::AbstractGeoTable, vars::AbstractVector)
  rb = @groupby(gt, :IJK)
  for v in vars
    rb = @combine(rb, {v} = mean({v}))
  end
  rb
end

reblock(df) = reblock(df, [Symbol(v) for v in names(vars) if !(v in (:IJK, "IJK"))])
reblock(df, ivar) = reblock(df, [Symbol(ivar)])

function downscale(geotable::AbstractGeoTable, factors; blks=true, ds=nothing)
  all(factors .== 1) && (return geotable)
  ds = isnothing(ds) ? get_spacing(geotable) : ds
  hs = ustrip.([d / 2 for d in ds])
  subd = ustrip.([d / x for (d, x) in zip(ds, factors)])
  isub = ustrip.([collect(y) .* subd for y in Iterators.product([0:(x - 1) for x in factors]...)])
  outfun = blks ? new_subblocks : new_subpts
  chunks = collect(partition(geotable, UniformPartition(Threads.nthreads(), false)))
  new_geoms = foldxt(vcat, Transducers.Map(x -> outfun(x, hs, isub, subd)), withprogress(chunks; interval=1e-3))

  new_geoms
end

function downscale(table::AbstractDataFrame, factors; blks=true, coords=[:x, :y, :z])
  all(factors .== 1) && (return table)
  ds = get_spacing(table)
  geotable = georef(table, coords)
  out = downscale(geotable, factors; blks=false, ds=ds)
  dom = domain(out)

  cpairs = mapreduce(vcat, enumerate(coords)) do (ic, c)
    cvals = [ustrip(to(centroid(dom, i))[ic]) for i in 1:nelements(dom)]
    c => cvals
  end

  hcat(DataFrame((; cpairs...)), DataFrame(values(out)))
end

function new_subblocks(geotable, hs, isub, subd)
  new_geoms = mapreduce(vcat, zip(geotable.geometry, geotable.IJK)) do (g, ijk)
    oi = ustrip.(to(centroid(g))) - hs
    origins = [oi + vi for vi in isub]
    subblks = mapreduce(vcat, origins) do o
      minmax = [(o[x], o[x] + subd[x]) for x in 1:length(subd)]
      pts = [Point(x...) for x in Iterators.product(minmax...)]
      Hexahedron(pts...)
    end
    georef((IJK=[ijk for x in eachindex(origins)],), subblks)
  end
  new_geoms
end

function new_subpts(geotable, hs, isub, subd)
  new_geoms = mapreduce(vcat, zip(geotable.geometry, geotable.IJK)) do (g, ijk)
    oi = ustrip.(to(centroid(g))) - hs
    origins = [oi + vi for vi in isub]
    subpts = mapreduce(vcat, origins) do o
      coords = [mean([o[x], o[x] + subd[x]]) for x in 1:length(subd)]
      Point(coords...)
    end
    georef((IJK=[ijk for x in eachindex(origins)],), subpts)
  end
  new_geoms
end

function to_dataframe(df)
  cols = names(df)
  rows = 1:nrow(df)
  sblks = DataFrame(values(df))
  if !("x" in cols)
    sblks[!, :x] = [ustrip(to(centroid(df.geometry, i))[1]) for i in rows]
    sblks[!, :y] = [ustrip(to(centroid(df.geometry, i))[2]) for i in rows]
    sblks[!, :z] = [ustrip(to(centroid(df.geometry, i))[3]) for i in rows]
  end
  sblks
end

function merge_subblocks(model1, model2)
  sblk1 = model1 isa DataFrame ? model1 : to_dataframe(model1)
  sblk1 = sblk1 |> Sort([:IJK, :x, :y, :z])

  sblk2 = model2 isa DataFrame ? model2 : to_dataframe(model2)
  sblk2 = sblk2 |> Sort([:IJK, :x, :y, :z])

  ijk1 = countmap(sblk1[:, :IJK])
  sub1 = [key for (key, val) in ijk1 if val > 1]
  reg1 = [key for (key, val) in ijk1 if val == 1]

  ijk2 = countmap(sblk2[!, :IJK])
  sub2 = [key for (key, val) in ijk2 if val > 1]
  reg2 = [key for (key, val) in ijk2 if val == 1]

  # this assumes no repeated columns except IJK; otherswise it may crash
  to_merge = []
  f1 = sub1 ∩ sub2
  println("- Merging sub-sub...")
  if length(f1) > 0
    i1 = [x in f1 for x in sblk1.IJK]
    i2 = [x in f1 for x in sblk2.IJK]
    push!(to_merge, hcat(sblk1[i1, :], sblk2[i2, Not(:IJK, :x, :y, :z)]))
  end

  f2 = reg1 ∩ reg2
  println("- Merging reg-reg...")
  if length(f2) > 0
    i1 = [x in f2 for x in sblk1.IJK]
    i2 = [x in f2 for x in sblk2.IJK]
    push!(to_merge, hcat(sblk1[i1, :], sblk2[i2, Not(:IJK, :x, :y, :z)]))
  end

  f3a = setdiff(sub1, f1)
  println("- Merging sub1-reg2...")
  if length(f3a) > 0
    i1 = [x in f3a for x in sblk1.IJK]
    i2 = [x in f3a for x in sblk2.IJK]
    newcols = setdiff(names(sblk2), names(sblk1))
    push!(to_merge, leftjoin(sblk1[i1, :], sblk2[i2, Not(:x, :y, :z)], on=:IJK))
  end

  f4a = setdiff(sub2, f1)
  println("- Merging sub2-reg1...")
  if length(f4a) > 0
    i1 = [x in f4a for x in sblk2.IJK]
    i2 = [x in f4a for x in sblk1.IJK]
    newcols = setdiff(names(sblk1), names(sblk2))
    push!(to_merge, leftjoin(sblk2[i1, :], sblk1[i2, Not(:x, :y, :z)], on=:IJK))
  end

  println("- Final vcat...")
  # merge final results
  sblk1 = vcat(to_merge...) |> Sort([:IJK, :x, :y, :z])

  sblk1
end

merge_subblocks(models...) = reduce(merge_subblocks, models)

# estimation case
function regblks_to_subblks(ptable, model, doms_filter; mode=:standard)
  suf = (mode in [:idw, :IDW]) ? "_IDW" : (mode in [:nn, :NN]) ? "_NN" : ""

  est_group = groupby(ptable, [:var, :model_domain])
  model[!, :idx] = 1:nrow(model)
  println("\nBLOCK TO SUBBLOCKS$suf")
  for ev in unique(ptable.var)
    ev == "RCATG" && mode in (:idw, :IDW, :nn, :NN) && continue
    model[!, "$ev$suf"] .= NaN
  end

  for group in est_group
    evar = group[!, :var][1]
    dom = group[!, :name_domain][1]
    evar == "RCATG" && mode in (:idw, :IDW, :nn, :NN) && continue
    println(" - Merging $evar $dom")
    emodel = CSV.read("regest_$(evar)_$(dom)$(suf).csv", DataFrame)

    filter_dom = doms_filter[group[!, :model_domain][1]]
    fmodel = filter(filter_dom, model)[!, ["idx", "IJK"]]
    fmodel = leftjoin(fmodel, emodel, on=:IJK) |> DropMissing()
    model[fmodel.idx, "$evar$suf"] .= fmodel[!, evar]
  end

  model[!, Not(:idx)]
end

# simulation case
function regblks_to_subblks(ptable, model, doms_filter, output)
  est_group = groupby(ptable, [:var, :model_domain])
  model[!, :idx] = 1:nrow(model)
  suffixes =
    output isa AbstractVector ? ["_q$q" for q in output] :
    output == "full" ? (1:maximum(ptable.nreals)) :
    output == "mean" ? ["_mean"] : error("output must be \"full\", \"mean\" or a vector of quantiles ranging [0,1]")

  println("\nBLOCK TO SUBBLOCKS")
  for ev in unique(ptable.var)
    for s in suffixes
      model[!, "$ev$s"] .= NaN
    end
  end

  for group in est_group
    evar = group[!, :var][1]
    dom = group[!, :name_domain][1]
    println(" - Merging $evar $dom")
    emodel = CSV.read("regsim_$(evar)_$(dom).csv", DataFrame)
    vars = [x for x in names(emodel) if x != "IJK"]
    if output isa AbstractVector
      colnames = ["$evar$s" for s in suffixes]
      for i in 1:length(output)
        emodel[!, colnames[i]] = [quantile(row, output[i]) for row in eachrow(emodel[!, vars])]
      end
      select!(emodel, "IJK", colnames...)
    elseif output == "mean"
      colnames = ["$(evar)_mean"]
      emodel[!, colnames[1]] = [mean(row) for row in eachrow(emodel[!, vars])]
      select!(emodel, "IJK", colnames...)
    else
      colnames = ["$evar$s" for s in suffixes]
    end

    filter_dom = doms_filter[group[!, :model_domain][1]]
    fmodel = filter(filter_dom, model)[!, ["idx", "IJK"]]
    fmodel = leftjoin(fmodel, emodel, on=:IJK) |> DropMissing()
    model[fmodel.idx, colnames] .= fmodel[!, colnames]
  end

  model[!, Not(:idx)]
end

function backflag(model, pts, varn, radius, blksize, sblksize)
  ijks = countmap(model.IJK)
  #sub = [key for (key, val) in ijks if val > 1]
  reg = [key for (key, val) in ijks if val == 1]

  ball = BallSearch(model.geometry, MetricBall(radius))

  inds = mapreduce(vcat, pts.geometry) do pt
    p = centroid(pt)
    search(p, ball)
  end

  inds = sort(unique(inds))
  backflag = typeof(getproperty(model, varn))(undef, nrow(pts))

  for blk in Tables.rows(model[inds, :])
    bsize = blk.IJK in reg ? blksize : sblksize
    o = ustrip.(to(centroid(blk.geometry))) - (bsize ./ 2)
    m = o + bsize
    bbox = Box(Point(o...), Point(m...))
    bf = findall(h -> centroid(h) ∈ bbox, pts.geometry)
    backflag[bf] .= blk[varn]
  end
  backflag
end
