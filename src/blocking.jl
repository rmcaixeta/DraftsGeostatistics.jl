
using GeoStats
using Transducers
using DataFrames
using StaticArrays

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

function backflag(model, pts, varn, radius, bsize)
  ball = BallSearch(model.geometry, MetricBall(radius))

  inds = mapreduce(vcat, pts.geometry) do pt
    p = centroid(pt)
    search(p, ball)
  end

  inds = sort(unique(inds))
  #backflag = typeof(getproperty(model, varn))(undef, nrow(pts))
  backflag = Vector{Float64}(undef, nrow(pts))

  for blk in Tables.rows(model[inds, :])
    o = ustrip.(to(centroid(blk.geometry))) - (bsize ./ 2)
    m = o + bsize
    bbox = Box(Point(o...), Point(m...))
    bf = findall(h -> centroid(h) ∈ bbox, pts.geometry)
    backflag[bf] .= blk[varn]
  end
  backflag
end

### refactored
function make_subblocks(
  centroids::AbstractMatrix,
  parentsize::Tuple,
  discretization::Tuple=(4, 4, 4);
  ijk=nothing,
  gtab=true
)
  nx, ny, nz = discretization
  n_sub = nx * ny * nz
  n_blocks = size(centroids, 1)
  ijks = isnothing(ijk) ? repeat(1:n_blocks, inner=n_sub) : repeat(ijk, inner=n_sub)

  x_offs = range(-0.5 + 1/(2*nx), 0.5 - 1/(2*nx), length=nx)
  y_offs = range(-0.5 + 1/(2*ny), 0.5 - 1/(2*ny), length=ny)
  z_offs = range(-0.5 + 1/(2*nz), 0.5 - 1/(2*nz), length=nz)

  xx = repeat(x_offs, inner=ny*nz)
  yy = repeat(y_offs, inner=nz, outer=nx)
  zz = repeat(z_offs, outer=nx*ny)
  offs = hcat(xx, yy, zz)
  out = reshape(centroids, n_blocks, 1, 3) .+ reshape(offs, 1, n_sub, 3) .* reshape(collect(parentsize), 1, 1, 3)

  cnames = [:x, :y, :z]
  out = (; zip(cnames, eachcol(reshape(out, :, 3)))..., :IJK => ijks)
  gtab ? georef(out, cnames) : out
end

function make_subblocks(domain::Domain, discretization::Tuple=(4, 4, 4); ijk=nothing)
  spacing = domain isa SubDomain ? parent(domain).spacing : domain.spacing
  parentsize = ustrip.(spacing)
  centroids = mapreduce(vcat, domain) do blk
    hcat(ustrip.(to(centroid(blk)))...)
  end
  make_subblocks(centroids, parentsize, discretization; ijk)
end

function make_subblocks(tab::AbstractGeoTable, discretization::Tuple=(4, 4, 4); ijk=:IJK)
  ijks = hasproperty(tab, ijk) ? getproperty(tab, ijk) : nothing
  make_subblocks(domain(tab), discretization; ijk=ijks)
end

# coord_to_ijk
# cartid(coord, c) = Int((coord - morigin[c]) ./ spacing_[c] .+ 1)
# ijks = LinearIndices((1:dims20[1],1:dims20[2],1:dims20[3]))
# ijks20c = [CartesianIndex(cartid(ref[i,:x],1), cartid(ref[i,:y],2), cartid(ref[i,:z],3)) for i in 1:nrow(ref)]
# ijks20 = [ijks[i] for i in ijks20c]
# ref.IJK = ijks20
