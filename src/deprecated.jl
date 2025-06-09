
function estim(group, evar, to_estim, allcomps, comps_filter, discr)
  if evar == :RCATG
    class_estim(group, evar, to_estim, allcomps, comps_filter)
  else
    grades_estim(group, evar, to_estim, allcomps, comps_filter, discr)
  end
end

function grades_estim(group, evar, to_estim, allcomps, comps_filter, discr)
  est = mapreduce(geovcat, eachrow(group)) do ln
    if to_estim == nothing
      to_estim
    else
      comps = allcomps |> GeoStats.Filter(comps_filter[parentindices(ln)[1]]) |> DropMissing(evar)
      # might need to read the comps here if different files are used
      estimator =
        isnothing(ln.estimator_args) || ln.estimator_args == "" ? evalstr(ln.estimator)() :
        evalstr(ln.estimator)(evalstr(ln.estimator_args))
      # need to deal with multi args
      nhood =
        ismissing(ln.neighborhood) ? nothing :
        startswith(ln.neighborhood, "AdvBallSearch") ?
        AdvBallSearch(comps; evalstr(replace(ln.neighborhood, "AdvBallSearch" => ""))...) : evalstr(ln.neighborhood)
      minneigh = eval(ln.minneighbors)
      maxneigh = eval(ln.maxneighbors)
      blkestim = ln.estimator == "Kriging" ? (point=false,) : (point=true,)

      discr = isnothing(discr) ? ntuple(i -> 4, embeddim(domain(to_estim))) : discr
      to_estim_ = ln.estimator == "Kriging" ? to_estim : downscale(to_estim, discr, blks=false)

      nblks = nrow(to_estim_)
      println("  - pass $(ln.pass), $nblks blks to estimate")

      partitions = collect(partition(to_estim_.geometry, UniformPartition(minimum([nthr, nblks]), false)))
      estfun = mapxt(
        gridchunk ->
          comps |>
          Select(evar) |>
          InterpolateNeighbors(
            gridchunk,
            model=estimator,
            minneighbors=minneigh,
            maxneighbors=maxneigh,
            neighborhood=nhood;
            blkestim...
          )
      )
      preest = foldxt(vcat, estfun, partitions)
      preest = hcat(preest, to_estim_)
      preest = ln.estimator == "Kriging" ? preest : georef(values(reblock(preest, [evar])), to_estim.geometry)

      missing_vals = findall(ismissing, getproperty(preest, evar))
      if length(missing_vals) > 0
        missing_ijks = to_estim.IJK[missing_vals]
        to_estim = to_estim |> GeoStats.Filter(row -> row.IJK in missing_ijks)
      else
        to_estim = nothing
        println("  - finished!")
      end
      preest |> DropMissing(evar)
    end
  end
end

function regblks_estimation(
  ptable,
  model,
  grid,
  allcomps,
  doms_filter,
  comps_filter;
  mode=:standard,
  idw=2,
  discr=nothing
)
  if mode in (:idw, :IDW)
    reftable = copy(ptable)
    reftable[!, :estimator] .= "IDW"
    reftable[!, :estimator_args] .= "$idw"
    suf = "_IDW"
  elseif mode in (:nn, :NN)
    reftable = copy(ptable)
    reftable[!, :estimator] .= "NN"
    reftable[!, :estimator_args] .= ""
    reftable[!, :neighborhood] .= ""
    reftable[!, :minneighbors] .= 1
    reftable[!, :maxneighbors] .= 1
    suf = "_NN"
  else
    reftable = ptable
    suf = ""
  end

  est_group = groupby(reftable, [:var, :model_domain])
  println("\nBLOCK ESTIMATES$suf")
  for group in est_group
    evar, dom = group[!, :var][1], group[!, :name_domain][1]
    evar == "RCATG" && mode in (:idw, :IDW, :nn, :NN) && continue

    filter_dom = doms_filter[group[!, :model_domain][1]]
    valid_ijks = unique(filter(filter_dom, model).IJK)

    to_estim = grid |> GeoStats.Filter(row -> row.IJK in valid_ijks)
    println("- $evar, $dom")
    est = estim(group, Symbol(evar), to_estim, allcomps, comps_filter, discr)
    CSV.write("regest_$(evar)_$(dom)$(suf).csv", values(est |> Sort(:IJK)))
  end
end

function class_estim(group, evar, to_estim, allcomps, comps_filter)
  est = mapreduce(geovcat, eachrow(group)) do ln
    if to_estim == nothing
      to_estim
    else
      nblks = nrow(to_estim)
      println("  - class $(ln.pass), $nblks blks to evaluate")
      comps = allcomps |> GeoStats.Filter(comps_filter[parentindices(ln)[1]])
      # might need to read the comps here if different files are used

      nhood =
        ismissing(ln.neighborhood) ? nothing :
        startswith(ln.neighborhood, "AdvBallSearch") ?
        AdvBallSearch(comps; evalstr(replace(ln.neighborhood, "AdvBallSearch" => ""))...) : evalstr(ln.neighborhood)
      minneigh = eval(ln.minneighbors)
      maxneigh = eval(ln.maxneighbors)
      # need to really consider or not the other parameters in here; currently only the neighborhood is used
      ctg = eval(ln.pass)

      partitions = collect(partition(to_estim, UniformPartition(minimum([nthr, nblks]), true)))
      estfun = mapxt(gridchunk -> class_blk(domain(comps), gridchunk, ctg, nhood))
      preest = foldxt(vcat, estfun, partitions)
      missing_vals = nrow(preest)
      preest = preest |> DropNaN(evar)
      missing_vals = missing_vals - nrow(preest)

      if missing_vals > 0
        to_estim = to_estim |> GeoStats.Filter(row -> !(row.IJK in preest.IJK))
      else
        to_estim = nothing
        println("  - finished!")
      end
      preest
    end
  end
end

function global_mean_validation(
  model::AbstractDataFrame,
  ptable::AbstractDataFrame,
  doms_filter::Dict,
  scenarios::AbstractVector,
  weightvar::Symbol
)
  name_domains = unique(ptable[!, [:model_domain, :name_domain]])
  vars = [x for x in unique(ptable.var) if x != "RCATG"]

  mean_checks = mapreduce(vcat, eachrow(name_domains)) do row
    dom = row.model_domain
    ndom = row.name_domain
    filter_dom = doms_filter[dom]
    out = filter(filter_dom, model)
    pairs = []
    for v in vars
      out[!, v] .*= out[!, weightvar]
      push!(pairs, v => sum => v)
      for s in scenarios
        out[!, "$(v)_$s"] .*= out[!, weightvar]
        push!(pairs, "$(v)_$s" => sum => "$(v)_$s")
      end
    end
    push!(pairs, weightvar => sum => weightvar)

    out = combine(groupby(out, :RCATG), pairs...)

    for v in vars
      out[!, v] ./= out[!, weightvar]
      for s in scenarios
        out[!, "$(v)_$s"] ./= out[!, weightvar]
        out[!, "$(v)_$(s)/estimate"] .= out[!, "$(v)_$s"] ./ out[!, v]
      end
    end
    insertcols!(out, 1, :domain => ndom)
    out |> Sort(:domain, :RCATG)
  end

  mean_checks
end

################################
############  SIMS  ############
################################

function simul(group, evar, to_estim, allcomps, comps_filter)
  grades_simul(group, evar, to_estim, allcomps, comps_filter)
end

function grades_simul(group, evar, to_estim, allcomps, comps_filter)
  ln = group[1, :]
  if isnothing(to_estim)
    nothing, nothing, nothing, nothing
  else
    nblks = nrow(to_estim)
    comps = allcomps |> GeoStats.Filter(comps_filter[parentindices(ln)[1]]) |> DropMissing(evar)
    wgts = ismissing(ln.wgt) ? nothing : weight(comps, evalstr(ln.wgt))
    hd_ns, refd = nscore(comps, evar; weights=wgts)

    # might need to read the comps here if different files are used
    vario = evalstr(ln.nsvariogram)
    nhood =
      ismissing(ln.neighborhood) ? nothing :
      startswith(ln.neighborhood, "AdvBallSearch") ?
      AdvBallSearch(comps; evalstr(replace(ln.neighborhood, "AdvBallSearch" => ""))...) : evalstr(ln.neighborhood)
    maxneigh = eval(ln.maxneighbors)
    nreals = eval(ln.nreals)
    nssill = eval(ln.nssill)

    vx = (nssill / sill(vario)) * vario
    process = GaussianProcess(vx)
    method = SEQMethod(neighborhood=nhood, maxneighbors=maxneigh)

    ## single
    sims = rand(StableRNG(2017), process, to_estim.geometry, hd_ns, nreals, method)

    ## multithread test
    #partitions = nreal_partitions(nreals, nthr)
    #simfun = mapxt(reals -> rand(StableRNG(2017*reals[1]), process, to_estim.geometry, hd_ns, length(reals), method))
    #sims = foldxt(simsvcat, simfun, partitions)

    ## multiprocessors test
    #sims = rand(StableRNG(2017), process, to_estim.geometry, hd_ns, nreals, method, pool=workers())

    # getting vario pts coords
    npts = nrow(to_estim)
    nsamps = 0.1 * npts < 10000 ? min(10000, npts) : min(round(Int, 0.1 * npts), 80000)
    idxs_for_vario = nsamps == npts ? (1:npts) : sample(1:npts, nsamps, replace=false)
    variopts = mapreduce(vcat, view(to_estim, idxs_for_vario).geometry) do geom
      coords = ustrip.(to(centroid(geom)))
      hcat(collect(coords)...)
    end
    coordnames = [:x, :y, :z]
    variopts = DataFrame(variopts, [coordnames[i] for i in 1:size(variopts, 2)])

    # getting histograms
    hist_qs = LinRange(0, 1, nrow(comps))
    histdf = DataFrame(q=hist_qs, samples=refd.values)

    mstats = sims_stats(sims, evar)
    println("  - avg of means: $(mstats[1,:mean]) | avg of vars: $(mstats[1,:var]) ...")

    println("  - reblocking ...")
    sims = mapreduce(hcat, 1:nreals) do i
      simvar = "$evar$i"
      #mx, vy = mean(sims[i].pct_Ni), var(sims[i].pct_Ni)
      sim = back_nscore(sims, i, evar, refd)
      variopts[!, simvar] .= sim[idxs_for_vario, evar]
      histdf[!, simvar] .= quantile(sim[!, evar], hist_qs)
      sim[!, :IJK] = to_estim.IJK
      sim = reblock(sim, [evar]) |> Rename(evar => simvar)
      i == 1 ? sim : sim[!, Not(:IJK)]
    end
    sims, histdf, mstats, variopts
  end
end

function sims_stats(sims::Ensemble, evar)
  reals = sims.reals
  mstats = mapreduce(vcat, reals) do sim
    hcat(mean(sim[evar]), var(sim[evar]))
  end
  mstats = round.(mean(mstats, dims=1), digits=2)
  DataFrame(mstats, [:mean, :var])
end

function sims_stats(sims::AbstractDataFrame)
  mstats = mapreduce(vcat, names(sims)) do col
    hcat(mean(sims[!, col]), var(sims[!, col]))
  end
  mstats = round.(mean(mstats, dims=1), digits=2)
  DataFrame(mstats, [:mean, :var])
end

function regblks_simulation(
  ptable,
  model,
  grid,
  allcomps,
  doms_filter,
  comps_filter,
  discr;
  histfolder="./sims_hist/",
  variofolder="./sims_vario/"
)
  est_group = groupby(ptable, [:var, :model_domain])
  println("\nBLOCK SIMULS")
  for group in est_group
    evar, dom = group[!, :var][1], group[!, :name_domain][1]
    println("- $evar, $dom")

    filter_dom = doms_filter[group[!, :model_domain][1]]
    valid_ijks = unique(filter(filter_dom, model).IJK)

    println("  - downscaling...")
    to_estim = grid |> GeoStats.Filter(row -> row.IJK in valid_ijks)
    to_estim = downscale(to_estim, discr, blks=false)

    println("  - simulating...")
    sims, histdf, mstats, variopts = simul(group, Symbol(evar), to_estim, allcomps, comps_filter)

    CSV.write("regsim_$(evar)_$(dom).csv", sims)
    CSV.write("$(histfolder)valsim_hist_$(evar)_$(dom).csv", histdf)
    CSV.write("$(histfolder)valsim_shist_$(evar)_$(dom).csv", mstats)
    CSV.write("$(variofolder)valsim_variog_$(evar)_$(dom).csv", variopts)
  end
end

# from declus
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
