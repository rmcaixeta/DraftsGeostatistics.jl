using GeoStats
using CSV
using DataFrames
using Transducers


nthr = Base.Threads.nthreads()
mapxt = Transducers.Map
geovcat(items...) = vcat([x for x in items if x != nothing]...)

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

            comps =
                allcomps |>
                GeoStats.Filter(comps_filter[parentindices(ln)[1]]) |>
                DropMissing(evar)
            # might need to read the comps here if different files are used
            estimator =
                isnothing(ln.estimator_args) || ln.estimator_args == "" ?
                evalstr(ln.estimator)() :
                evalstr(ln.estimator)(evalstr(ln.estimator_args))
            # need to deal with multi args
            nhood =
                ismissing(ln.neighborhood) ? nothing :
                startswith(ln.neighborhood, "AdvBallSearch") ?
                AdvBallSearch(
                    comps;
                    evalstr(replace(ln.neighborhood, "AdvBallSearch" => ""))...,
                ) : evalstr(ln.neighborhood)
            minneigh = eval(ln.minneighbors)
            maxneigh = eval(ln.maxneighbors)
            blkestim = ln.estimator == "Kriging" ? (point = false,) : (point = true,)

            discr =
                isnothing(discr) ? ntuple(i -> 4, embeddim(domain(to_estim))) : discr
            to_estim_ =
                ln.estimator == "Kriging" ? to_estim :
                downscale(to_estim, discr, blks = false)

            nblks = nrow(to_estim_)
            println("  - pass $(ln.pass), $nblks blks to estimate")

            partitions = collect(
                partition(
                    to_estim_.geometry,
                    UniformPartition(minimum([nthr, nblks]), false),
                ),
            )
            estfun = mapxt(
                gridchunk ->
                    comps |> InterpolateNeighbors(
                        gridchunk,
                        evar => estimator,
                        minneighbors = minneigh,
                        maxneighbors = maxneigh,
                        neighborhood = nhood;
                        blkestim...,
                    ),
            )
            preest = foldxt(vcat, estfun, partitions)
            preest = hcat(preest, to_estim_)
            preest =
                ln.estimator == "Kriging" ? preest :
                georef(values(reblock(preest, [evar])), to_estim.geometry)

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
    mode = :standard,
    idw = 2,
    discr = nothing,
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
                AdvBallSearch(
                    comps;
                    evalstr(replace(ln.neighborhood, "AdvBallSearch" => ""))...,
                ) : evalstr(ln.neighborhood)
            minneigh = eval(ln.minneighbors)
            maxneigh = eval(ln.maxneighbors)
            # need to really consider or not the other parameters in here; currently only the neighborhood is used
            ctg = eval(ln.pass)

            partitions = collect(
                partition(to_estim, UniformPartition(minimum([nthr, nblks]), true)),
            )
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



function class_blk(dh, to_estim, ct, advsearch)
    grid = to_estim.geometry
    knn = KNearestSearch(dh, 1)
    rangedist = ustrip(advsearch.ball.radii[1])
    catgs = mapreduce(vcat, grid) do blk
        p = centroid(blk)
        n = search(p, knn)
        d = evaluate(Euclidean(), ustrip.(to(p)), ustrip.(to(centroid(dh, n[1]))))
        outval = NaN

        if d <= rangedist
            n = search(p, advsearch)
            if length(n) == 3
                dists = [
                    evaluate(Euclidean(), ustrip.(to(p)), ustrip.(to(centroid(dh, ind)))) for ind in n
                ]
                outval = (sum(dists) < rangedist * 2) ? ct : NaN
            end
        end
        outval
    end

    georef((RCATG = catgs, IJK = to_estim.IJK), grid)
end



function global_mean_validation(
    model::AbstractDataFrame,
    ptable::AbstractDataFrame,
    doms_filter::Dict,
    scenarios::AbstractVector,
    weightvar::Symbol,
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
