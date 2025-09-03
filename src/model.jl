
## MAIN WORKFLOW: categorical_compositing() -> extract_intrusion_pts() -> intrusion_model()

function categorical_compositing(
  dh_,
  var,
  interior,
  exterior;
  code_var=:code_catg_,
  info_var=:info_catg_,
  filters=[("Exterior", 5.0), ("Interior", 2.1)]
)
  pars = dh_.pars
  bh_, fr_, to_ = pars.holeid, pars.from, pars.to
  dh = dh_.table |> Sort(bh_, fr_)

  # filters
  nvals = nrow(dh)
  interior_filter = findall(x -> x in interior, dh[!, var])
  exterior_filter = findall(x -> x in exterior, dh[!, var])

  catgtype = ["Ignored" for i in 1:nvals]
  catgtype[interior_filter] .= "Interior"
  catgtype[exterior_filter] .= "Exterior"

  # table grouped by codes
  code = set_interval_code_(dh, catgtype, bh_)
  dht = hcat(dh, DataFrame((; code_var=>code, info_var=>catgtype)))
  grp = groupby(dht, [code_var, info_var, bh_])
  grp = combine(grp, fr_ => minimum => fr_, to_ => maximum => to_)
  grp.LENGTH = grp[!, to_] - grp[!, fr_]

  # return here if ignored is wanted and no filters of intervals

  ## ignore ignored
  grp = grp[grp[!, info_var] .!= "Ignored", :]

  # get upper and lower bounds of the hole, to not filter
  holeids = unique(grp[!, bh_])
  bounds =
    vcat([searchsortedfirst(grp[!, bh_], x) for x in holeids], [searchsortedlast(grp[!, bh_], x) for x in holeids])

  # apply filters for small intervals
  for (ctype, cleng) in filters
    filter_t = (grp[!, info_var] .== ctype) .& (grp[!, :LENGTH] .< cleng)
    filter_t = [x for x in findall(filter_t) if !(x in bounds)]

    filter_t = [x for x in 1:nrow(grp) if !(x in filter_t)]
    grp = grp[filter_t, :]
    code = set_interval_code_(grp, grp[!, info_var], bh_)
    grp = hcat(select(grp, Not(code_var)), DataFrame((; code_var=>code)))

    grp = groupby(grp, [code_var, info_var, bh_])
    grp = combine(grp, fr_ => minimum => fr_, to_ => maximum => to_)
    grp.LENGTH = grp[!, to_] - grp[!, fr_]
  end
  fillxyz!(grp, dh_.trace, pars, output=["mid", "from", "to"])
  DrillHole(grp, dh_.trace, pars, dh_.warns)
end

function set_interval_code_(dh, catg, bh_)
  # aux function; this loop does not consider gaps currently
  # dh must be already sorted
  code = [1]
  for i in 2:nrow(dh)
    c = code[end]
    if dh[i - 1, bh_] == dh[i, bh_] && catg[i - 1] == catg[i]
      push!(code, c)
    else
      push!(code, c + 1)
    end
  end
  code
end

function signed_distances(gtab::AbstractGeoTable, var=:info_catg_, interior=["Interior"], exterior=["Exterior"])
  codes = getproperty(gtab, Symbol(var))
  int_filter = findall(x -> x in interior, codes)
  ext_filter = findall(x -> x in exterior, codes)
  dom = domain(gtab)

  s = KNearestSearch(view(dom, ext_filter), 1)
  interior_sd = mapreduce(vcat, int_filter) do i
    _, d = searchdists(centroid(dom, i), s)
    -ustrip(d[0])
  end

  s = KNearestSearch(view(dom, int_filter), 1)
  exterior_sd = mapreduce(vcat, ext_filter) do i
    _, d = searchdists(centroid(dom, i), s)
    ustrip(d[0])
  end

  sd = fill(NaN, length(dom))
  sd[int_filter] .= interior_sd
  sd[ext_filter] .= exterior_sd
  gsd = georef((; :SD => sd), dom)
  hcat(gtab, gsd)
end

function extract_intrusion_pts(dh_, var=:info_catg_, interior=["Interior"], exterior=["Exterior"]; spacing=1.0)
  dh = copy(dh_.table)
  bh_, fr_ = dh_.pars.holeid, dh_.pars.from

  interior_filter = findall(x -> x in interior, dh[!, var])
  exterior_filter = findall(x -> x in exterior, dh[!, var])
  bhids = dh[!, bh_]

  # split intervals in spacing m points or minimum of two values per interval; this rule can be improved
  interior_cmp = mapreduce(vcat, interior_filter) do i
    lngt, bh = abs(dh[i, :LENGTH]), dh[i, bh_]
    fr = dh[i, fr_]
    parts = ceil(Int, lngt / spacing) + 1
    intervals = LinRange(0.01, lngt - 0.01, parts)
    ofrom = [fr + itx for itx in intervals]
    obhid = [bh for itx in intervals]
    tab = DataFrame(bh_ => obhid, fr_ => ofrom)
    fillxyz!(tab, dh_.trace, dh_.pars, output=["from"])
    [(row[bh_], row[:X_FROM], row[:Y_FROM], row[:Z_FROM]) for row in eachrow(tab)]
  end

  exterior_cmp = mapreduce(vcat, exterior_filter) do i
    lngt, bh = abs(dh[i, :LENGTH]), dh[i, bh_]
    fr = dh[i, fr_]
    parts = ceil(Int, lngt / spacing) + 1
    intervals = LinRange(0.01, lngt - 0.01, parts)
    ofrom = [fr + itx for itx in intervals]
    obhid = [bh for itx in intervals]
    tab = DataFrame(bh_ => obhid, fr_ => ofrom)
    fillxyz!(tab, dh_.trace, dh_.pars, output=["from"])
    [(row[bh_], row[:X_FROM], row[:Y_FROM], row[:Z_FROM]) for row in eachrow(tab)]
  end

  geoms = vcat([Point(x[2:end]...) for x in interior_cmp], [Point(x[2:end]...) for x in exterior_cmp])
  bhids = vcat([x[1] for x in interior_cmp], [x[1] for x in exterior_cmp])
  vals = vcat(["Interior" for x in interior_cmp], ["Exterior" for x in exterior_cmp])
  tab = (; bh_ => bhids, var => vals)
  georef(tab, PointSet(geoms))
end

LocalEstimator = Union{LocalKrigingModel,LocalIDWModel}

function intrusion_model(
  sd,
  interpolant,
  grid;
  subblocks_split=(3, 3, 8),
  sub_interpolant=IDW(3),
  maxneighbors=100,
  sub_maxneighbors=50
)
  sdev = std(sd.SD)
  dh = @transform(sd, :SD = :SD / sdev)
  res_init = get_spacing(grid)

  # do a 1st round estimate in upscaled blocks
  ##

  # estimate in parent blocks
  out = grid isa CartesianGrid ? grid : grid.geometry
  interpfun = interpolant isa LocalEstimator ? LocalInterpolate : InterpolateNeighbors
  est = dh |> interpfun(out, :SD => interpolant, maxneighbors=maxneighbors)
  #add ijk
  est = hcat(est, georef((IJK=1:nrow(est),), est.geometry))
  lim = ustrip(maximum(res_init) / (sdev * 2))
  #filter
  m = ustrip.(abs.(est.SD)) .< lim
  ok = est[.!m, :]

  # estimate in subblocks
  ##downscale
  refined = downscale(est[m, :], subblocks_split)
  ##interpolate2

  est = if sub_interpolant isa LocalEstimator
    lp_ = nnpars(sub_interpolant.localaniso, grid, refined)
    sub_interpolant_ =
      sub_interpolant isa LocalKrigingModel ? LocalKriging(sub_interpolant.method, lp_, sub_interpolant.γ) :
      LocalIDW(sub_interpolant.exponent, lp_)
    dh |> Select(:SD) |> LocalInterpolate(refined.geometry, model=sub_interpolant_, maxneighbors=sub_maxneighbors)
  else
    chunks = collect(partition(refined.geometry, UniformPartition(Threads.nthreads(), false)))
    foldxt(
      vcat,
      Transducers.Map(
        x -> dh |> Select(:SD) |> InterpolateNeighbors(x, model=sub_interpolant, maxneighbors=sub_maxneighbors)
      ),
      chunks
    )
  end

  est = hcat(est, georef((IJK=refined.IJK,), refined.geometry))
  est = vcat(est, ok)
  est = @transform(est, :SD = :SD * sdev)
  est
end

function extract_contacts(grp; valid_holeids=nothing, code_var=:code_catg_, info_var=:info_catg_)
  ## i guess it works only if the stratigraphy is exterior at the bottom and interior at the top
  pars = grp.pars
  bh_, fr_, to_ = pars.holeid, pars.from, pars.to
  maxz = maximum(vcat(grp.table[!, :Z_FROM], grp.table[!, :Z_TO])) + 1
  minz = minimum(vcat(grp.table[!, :Z_FROM], grp.table[!, :Z_TO])) - 1

  dh_contact = groupby(grp.table, [bh_, info_var])
  dh_contact = combine(
    dh_contact,
    fr_ => minimum => fr_,
    to_ => maximum => to_,
    :Z_FROM => first => :MinZ,
    :Z_TO => last => :MaxZ,
    :Z_FROM => maximum => :Z0
  )
  dh_contact.MinZ .-= 0.01
  dh_contact.MaxZ .+= 0.01

  ## The folowing might be more correct, but need tests
  #dh_contact = combine(dh_contact, fr_ => minimum => fr_, to_ => maximum => to_, :Z_FROM => minimum => :MinZ_FROM, :Z_FROM => maximum => :MaxZ_FROM, :Z_TO => minimum => :MinZ_TO, :Z_TO => maximum => :MaxZ_TO)
  #dh_contact.Z0 = maximum(df[!, [:MaxZ_FROM, :MaxZ_TO]], dims=2) #[:]
  #dh_contact.MinZ = maximum(df[!, [:MaxZ_FROM, :MaxZ_TO]], dims=2) .- 0.01
  #dh_contact.MaxZ = minimum(df[!, [:MinZ_FROM, :MinZ_TO]], dims=2) .+ 0.01

  interior = dh_contact[!, info_var] .== "Interior"
  dh_contact[interior, :MinZ] .= minz
  dh_contact[.!interior, :MaxZ] .= maxz
  dh_contact = groupby(dh_contact, bh_)
  dh_contact = combine(
    dh_contact,
    fr_ => minimum => fr_,
    to_ => maximum => to_,
    :MinZ => maximum => :MinZ,
    :MaxZ => minimum => :MaxZ,
    info_var => first => :First,
    info_var => last => :Last,
    info_var => length => :COUNT,
    :Z0 => maximum => :Z0
  )

  contact_rule = (dh_contact[!, :COUNT] .> 1) .& (dh_contact[!, :First] .!= dh_contact[!, :Last])
  if !isnothing(valid_holeids)
    valid_holeids = [x in valid_holeids for x in dh_contact[!, bh_]]
    contact_rule = contact_rule .& valid_holeids
  end
  holeids_contact = dh_contact[contact_rule, bh_]
  holeids_nocontact = dh_contact[.!contact_rule, bh_]

  # Get Z contact; add some rule later to get first contact or last contact
  contact_pts = filter(row -> row[bh_] in holeids_contact, grp.table)
  contact_pts = groupby(contact_pts, [bh_, :info_catg_])
  contact_pts = combine(contact_pts, code_var => last => :code)

  contact_pts = contact_pts[contact_pts[!, :info_catg_] .== "Exterior", :code]
  contact_pts = filter(row -> row[code_var] in contact_pts, grp.table)
  contact_pts = select(contact_pts, bh_, :X_FROM => :X, :Y_FROM => :Y, :Z_FROM => :Z)

  # Merge info
  dh_contact = leftjoin(dh_contact, contact_pts, on=bh_)
  notmiss = .!ismissing.(dh_contact.Z)
  dh_contact[notmiss, :MinZ] .= dh_contact[notmiss, :Z]
  dh_contact[notmiss, :MaxZ] .= dh_contact[notmiss, :Z]
  badrule = dh_contact.MinZ .> dh_contact.MaxZ
  dh_contact[badrule, :MinZ] .= minz
  dh_contact[badrule, :MaxZ] .= maxz
  select!(dh_contact, Not(:COUNT))
  dh_contact
end

function impute_missing_contacts(contact_table, interpolant, comps; only_2d=false)
  pars = comps.pars
  bh_, fr_, to_, dp_ = pars.holeid, pars.from, pars.to, pars.dip
  maxz = maximum(contact_table.MaxZ)

  no_contact = filter(row -> ismissing(row.Z), contact_table)[!, bh_]
  to_estim = combine(groupby(comps.trace, bh_), dp_ => mean => dp_, :X => mean => :X, :Y => mean => :Y)
  to_estim = to_estim[abs.(to_estim[!, dp_]) .== 90, [bh_, :X, :Y]]
  to_estim = filter(row -> row[bh_] in no_contact, to_estim)
  hd = filter(row -> !ismissing(row.Z), contact_table)

  (nrow(to_estim) == 0 || nrow(hd) == 0) && (return comps)

  to_estim = leftjoin(to_estim, select(contact_table, Not(:X, :Y, :Z)), on=bh_)
  to_estim = georef(to_estim, (:X, :Y))
  hd = georef(hd, (:X, :Y))

  est = hd |> UniqueCoords(:Z => mean) |> InterpolateNeighbors(to_estim.geometry, :Z => interpolant, maxneighbors=100)
  est = hcat(to_estim, est)
  only_2d && (return est)

  newrows = @chain est begin
    @transform(:Z = clamp(ustrip.((:Z, :MinZ, :MaxZ))...))
    @transform(:AT = :Z0 - :Z)
    @transform(
      :info_catg_ = ifelse(isequal(:MaxZ, maxz), "Exterior", "Interior"),
      {fr_} = ifelse(isequal(:MaxZ, maxz), :AT, {to_}),
      {to_} = ifelse(isequal(:MaxZ, maxz), {fr_}, :AT)
    )
    @transform(:LENGTH = {to_} - {fr_})
  end

  newrows = newrows |> Select(:info_catg_, bh_, fr_, to_, :LENGTH, :Z0)
  newrows2 = @chain newrows begin
    @transform(:info_catg_ = ifelse(:info_catg_ == "Interior", "Exterior", "Interior"))
    @transform({fr_} = ifelse(:info_catg_ == "Interior", {fr_} - 10, {to_}))
    @transform({to_} = {fr_} + 10, :LENGTH = 10)
  end

  newrows = vcat(newrows, newrows2) |> Sort([bh_, fr_])
  newrows = DataFrame(values(newrows))

  newdh = vcat(comps.table, newrows, cols=:intersect) |> Sort([bh_, fr_])
  newdh[!, :code_catg_] = set_interval_code_(newdh, newdh[!, :info_catg_], bh_)
  newdh = groupby(newdh, [:code_catg_, :info_catg_, bh_])
  newdh = combine(newdh, fr_ => minimum => fr_, to_ => maximum => to_)
  newdh[!, :LENGTH] = newdh[!, to_] - newdh[!, fr_]

  fillxyz!(newdh, comps.trace, pars, output=["mid", "from", "to"])
  DrillHole(newdh, comps.trace, pars, comps.warns)
end

function extrapolate_borders(
  dh_;
  min_abs_dip=45.0,
  ignore_lower=[],
  ignore_upper=[],
  code_var=:code_catg_,
  info_var=:info_catg_
)
  pars = dh_.pars
  dip_negative = !pars.invdip

  bh_, fr_, to_, dp_ = pars.holeid, pars.from, pars.to, pars.dip
  dh = dh_.table |> Sort([bh_, fr_])
  maxz = maximum(dh.Z_FROM)
  minz = minimum(dh.Z_TO)

  dip_negative && (dh_valid[!, dp_] .*= -1)
  dh_valid = combine(groupby(dh_.trace, bh_), dp_ => first => dp_)
  dh_valid = dh_valid[abs.(dh_valid[!, dp_]) .> min_abs_dip, bh_]
  dh_valid = filter(row -> row[bh_] in dh_valid && !(row[bh_] in ignore_upper), dh)

  upper_extrap = groupby(dh_valid, bh_)
  upper_extrap = combine(
    upper_extrap,
    code_var => first => code_var,
    info_var => first => info_var,
    fr_ => first => to_,
    :Z_FROM => first => :Z
  )
  upper_extrap[!, :LENGTH] = maxz .- upper_extrap[!, :Z]
  upper_extrap[!, fr_] = upper_extrap[!, to_] - upper_extrap[!, :LENGTH]
  upper_extrap = upper_extrap[upper_extrap[!, info_var] .== "Interior", Not(:Z)]

  dh_valid = combine(groupby(dh_.trace, bh_), dp_ => last => dp_)
  dh_valid = dh_valid[abs.(dh_valid[!, dp_]) .> min_abs_dip, bh_]
  dh_valid = filter(row -> row[bh_] in dh_valid && !(row[bh_] in ignore_lower), dh)

  lower_extrap = groupby(dh_valid, bh_)
  lower_extrap = combine(
    lower_extrap,
    code_var => last => code_var,
    info_var => last => info_var,
    to_ => last => fr_,
    :Z_TO => last => :Z
  )
  lower_extrap[!, :LENGTH] = lower_extrap[!, :Z] .- minz
  lower_extrap[!, to_] = lower_extrap[!, fr_] + lower_extrap[!, :LENGTH]
  lower_extrap = lower_extrap[lower_extrap[!, info_var] .== "Exterior", Not(:Z)]

  newdh = vcat(dh, upper_extrap, lower_extrap, cols=:intersect) |> Sort([bh_, fr_])
  newdh[!, code_var] = set_interval_code_(newdh, newdh[!, info_var], bh_)
  newdh = groupby(newdh, [code_var, info_var, bh_])
  newdh = combine(newdh, fr_ => minimum => fr_, to_ => maximum => to_)
  newdh[!, :LENGTH] = newdh[!, to_] - newdh[!, fr_]

  fillxyz!(newdh, dh_.trace, pars, output=["mid", "from", "to"])
  DrillHole(newdh, dh_.trace, pars, dh_.warns)
end

function localaniso_from_pts(comps::DrillHole; ratios=[1.0, 1.0, 0.5], ball=200)
  holeid = comps.pars.holeid
  f = row -> -1 <= row[:SD] < 0
  sd = extract_intrusion_pts(comps, :info_catg_, ["Interior"], ["Exterior"]) |> GeoStats.Filter(f)

  searcher = AdvBallSearch(
    sd,
    MetricBall(ball),
    k=3,
    usesectors=(max=1, n=8, split=false),
    maxpercategory=(; holeid => 1),
    rank_metric=:same
  )
  geom = sd.geometry

  tri = mapreduce(vcat, geom) do pt
    p = centroid(pt)
    n = search(p, searcher)
    length(n) == 3 ? Triangle(geom[n]...) : missing
  end

  tri = filter(!ismissing, tri)
  tri = GeometrySet(tri)
  tri, localanisotropies(Geometric, tri, ratios)
end

function localaniso_from_pts(pts::AbstractGeoTable, holeid; ratios=[1.0, 1.0, 0.5], ball=200)
  f = row -> -1 <= row[:SD] < 0
  sd = pts |> GeoStats.Filter(f)

  searcher = AdvBallSearch(
    sd,
    MetricBall(ball),
    k=3,
    usesectors=(max=1, n=8, split=false),
    maxpercategory=(; holeid => 1),
    rank_metric=:same
  )
  geom = sd.geometry

  tri = mapreduce(vcat, geom) do pt
    p = centroid(pt)
    n = search(p, searcher)
    length(n) == 3 ? Triangle(geom[n]...) : missing
  end

  tri = filter(!ismissing, tri)
  tri = GeometrySet(tri)
  tri, localanisotropies(Geometric, tri, ratios)
end
