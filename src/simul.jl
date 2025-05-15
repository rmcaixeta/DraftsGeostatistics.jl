
nthr = Base.Threads.nthreads()
mapxt = Transducers.Map


function nreal_partitions(nreal, npart)
    flat_v = 1:nreal
    base_size, remainder = divrem(nreal, npart)
    v = Vector{Vector{Int}}()

    start_idx = 1
    for i = 1:npart
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


function simul(group, evar, to_estim, allcomps, comps_filter)
    grades_simul(group, evar, to_estim, allcomps, comps_filter)
end

function grades_simul(group, evar, to_estim, allcomps, comps_filter)
    ln = group[1, :]
    if isnothing(to_estim)
        nothing, nothing, nothing, nothing
    else
        nblks = nrow(to_estim)
        comps =
            allcomps |>
            GeoStats.Filter(comps_filter[parentindices(ln)[1]]) |>
            DropMissing(evar)
        wgts = ismissing(ln.wgt) ? nothing : weight(comps, evalstr(ln.wgt))
        hd_ns, refd = nscore(comps, evar; weights = wgts)

        # might need to read the comps here if different files are used
        vario = evalstr(ln.nsvariogram)
        nhood =
            ismissing(ln.neighborhood) ? nothing :
            startswith(ln.neighborhood, "AdvBallSearch") ?
            AdvBallSearch(
                comps;
                evalstr(replace(ln.neighborhood, "AdvBallSearch" => ""))...,
            ) : evalstr(ln.neighborhood)
        maxneigh = eval(ln.maxneighbors)
        nreals = eval(ln.nreals)
        nssill = eval(ln.nssill)

        vx = (nssill / sill(vario)) * vario
        process = GaussianProcess(vx)
        method = SEQMethod(neighborhood = nhood, maxneighbors = maxneigh)

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
        idxs_for_vario = nsamps == npts ? (1:npts) : sample(1:npts, nsamps, replace = false)
        variopts = mapreduce(vcat, view(to_estim, idxs_for_vario).geometry) do geom
            coords = ustrip.(to(centroid(geom)))
            hcat(collect(coords)...)
        end
        coordnames = [:x, :y, :z]
        variopts = DataFrame(variopts, [coordnames[i] for i = 1:size(variopts, 2)])

        # getting histograms
        hist_qs = LinRange(0, 1, nrow(comps))
        histdf = DataFrame(q = hist_qs, samples = refd.values)

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
    mstats = round.(mean(mstats, dims = 1), digits = 2)
    DataFrame(mstats, [:mean, :var])
end

function sims_stats(sims::AbstractDataFrame)
    mstats = mapreduce(vcat, names(sims)) do col
        hcat(mean(sims[!, col]), var(sims[!, col]))
    end
    mstats = round.(mean(mstats, dims = 1), digits = 2)
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
    histfolder = "./sims_hist/",
    variofolder = "./sims_vario/",
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
        to_estim = downscale(to_estim, discr, blks = false)

        println("  - simulating...")
        sims, histdf, mstats, variopts =
            simul(group, Symbol(evar), to_estim, allcomps, comps_filter)

        CSV.write("regsim_$(evar)_$(dom).csv", sims)
        CSV.write("$(histfolder)valsim_hist_$(evar)_$(dom).csv", histdf)
        CSV.write("$(histfolder)valsim_shist_$(evar)_$(dom).csv", mstats)
        CSV.write("$(variofolder)valsim_variog_$(evar)_$(dom).csv", variopts)
    end
end


# wgts = weight(comps, BlockWeighting(55,55,14))
function nscore(gtab, evar; weights = nothing)
    data = getproperty(gtab, evar)
    weights = isnothing(weights) || weights isa GeoWeights ? weights : GeoWeights(domain(gtab), weights)
    outvec, refd = nscore(data; weights)
    outtab = georef((; evar => outvec), gtab.geometry)
    outtab, refd
end

function nscore(data; weights = nothing)
    so = TableTransforms.qsmooth(data)
    ss =
        isnothing(weights) ? so :
        TableTransforms.qsmooth(quantile(data, weights, LinRange(0, 1, length(data))))
    refd = TableTransforms.EmpiricalDistribution(ss)
    TableTransforms.qtransform(so, refd, Normal()), refd
end

function back_nscore(sims::Ensemble, ireal::Int, evar, refd; georef = false)
    vals = sims.reals[ireal][evar]
    georef ? georef((; evar => back_nscore(vals, refd)), sims.domain) :
    DataFrame((; evar => back_nscore(vals, refd)))
end

function back_nscore(gtab::AbstractGeoTable, evar, refd; as_geotable = true)
    vals = getproperty(gtab, evar)
    as_geotable ? georef((; evar => back_nscore(vals, refd)), domain(gtab)) :
    DataFrame((; evar => back_nscore(vals, refd)))
end

back_nscore(vals::AbstractVector, refd) = TableTransforms.qtransform(vals, Normal(), refd)


#function validation_graphs(ptable)
#	est_group = groupby(ptable, [:var,:model_domain])
#
#	for group in est_group
#		evar = group[!, :var][1]
#		dom  = group[!, :name_domain][1]
#		println(" - Validations $evar $dom")
#
#		# "./sims_hist/valsim_hist_$(evar)_$(dom).csv", histdf
#		# "./sims_hist/valsim_shist_$(evar)_$(dom).csv", mstats
#		# "./sims_vario/valsim_variog_$(evar)_$(dom).csv", variopts
#
#		hmodel = CSV.read("./sims_hist/valsim_hist_$(evar)_$(dom).csv", DataFrame)
#		nstats = CSV.read("./sims_hist/valsim_shist_$(evar)_$(dom).csv", DataFrame)
#		vmodel = CSV.read("./sims_vario/valsim_variog_$(evar)_$(dom).csv", DataFrame)
#
#		# hist graphs
#		nu, nv = nstats[1,:mean], nstats[1,:var]
#		simcols = [x for x in names(hmodel) if !(x in ["samples","q"])]
#		fig = Mke.Figure()
#		ax = Mke.Axis(fig[1,1],xlabel=evar,ylabel="cdf",title="$(evar) - $(dom) ($nu,$nv)")
#		for col in simcols
#			Mke.lines!(ax, hmodel[!,col], hmodel.q, color=(:gray,0.6))
#		end
#		Mke.lines!(ax, hmodel.samples, hmodel.q, color=:black, linestyle=:dash)
#		Mke.save("./sims_hist/hist_$(evar)_$(dom).png", fig)
#
#		# vario graphs
#		vario = evalstr(group[!, :variogram][1])
#		mstats = sims_stats(vmodel[:,simcols])
#		dsill = mstats[1,:var]
#		vx = (dsill/sill(vario)) * vario
#		emps = [(fun=PlanarVariogram, args=((0.,0.,1.),), kwargs=(ntol=2.0, maxlag=160, nlags=80)),(fun=DirectionalVariogram, args=((0.,0.,1.),), kwargs=(maxlag=24,nlags=12,dtol=1.0))]
#		dirs = mapreduce(vcat, emps) do exp
#			if exp.fun == PlanarVariogram
#				n = collect(exp.args[1])
#				D1 = length(n)-1
#				v = [0 for x in n]
#				v[1] = n[1] == 0 ? 1 : v[1]
#				v[end] = -dot(n[1:D1],v[1:D1])/n[end]
#				v1 = v/norm(v)
#				v2 = cross(v1,n)
#				Tuple(v1)
#			elseif exp.fun == DirectionalVariogram
#				exp.args[1]
#			end
#		end
#
#		vmodel = georef(vmodel, (:x,:y,:z))
#
#
#		simsvar = mapreduce(vcat,simcols) do col
#			mapreduce(vcat,enumerate(emps)) do (i,emp)
#				name = "$(col)$(i)"
#				println("Doing $(col) $(dirs[i]) vario ...")
#				args1 = haskey(emp, :args) ? emp.args : []
#				args2 = (vmodel, col) # modify for multivar case
#				evario = emp.fun(args1..., args2...; emp.kwargs...)
#				x, y, p = values(evario)
#				n = [name for k in x]
#				hcat(n,ustrip.(x),y,p)
#			end
#		end
#
#		println("Doing theoretical varios ...")
#		modelsvar = mapreduce(vcat,enumerate(emps)) do (i,emp)
#			vfrom = Point(0,0,0)
#			pto = dirs[i]
#			maxlag = emp.kwargs.maxlag
#			mapreduce(vcat,1:maxlag) do x
#				vto = pto .* x
#				y = vx(vfrom, Point(vto...))
#				hcat("model$i",x,y,0)
#			end
#		end
#
#		varios = DataFrame(vcat(simsvar,modelsvar),[:name,:x,:y,:p])
#
#		# Loop through each directory
#		for d in 1:length(dirs)
#			endval = string(d)
#			dfx = filter(row -> endswith(row.name, endval), varios)
#			fig = Mke.Figure()
#			ax = Mke.Axis(fig[1, 1], title ="$(evar) - $(dom) - $(dirs[d])", xlabel = "Lag", ylabel = "Variogram")
#			for catg in unique(dfx.name)
#				dfy = filter(row -> row.name == catg, dfx)
#				if startswith(catg, "model")
#					Mke.lines!(ax, dfy.x, dfy.y)
#				else
#					Mke.lines!(ax, dfy.x, dfy.y, color = :gray, linestyle = :dash, alpha = 0.5)  # Other data, dashed gray lines
#				end
#			end
#
#			# Add center line for non-model data
#			dfe = filter(row -> !startswith(row.name, "model"), dfx)
#			dfe = combine(groupby(dfe, :x), :y => median => :y_median)  # Pivot-like behavior
#			Mke.lines!(ax, dfe.x, dfe.y_median, color = :black, linestyle = :dash, alpha = 0.8)
#			Mke.save("./sims_vario/vario_$(evar)_$(dom)_$(endval).png", fig)  # Save the figure to a PNG file
#		end
#
#	end
#end
#
#
#
#
#
#import pandas as pd
#import matplotlib.pyplot as plt
#import numpy as np
#varios = pd.read_csv('exp_results.csv')
#dirs = 2
#
### need to add center line
#
#for d in range(2):
#	endval = f'{d+1}'
#	dfx = varios[varios.name.str.endswith(endval)]
#	for catg in dfx.name.unique():
#		dfy = dfx[dfx.name == catg]
#		if catg.startswith('model'): plt.plot(dfy.x,dfy.y)
#		else: plt.plot(dfy.x,dfy.y,color='gray',alpha=0.5,linestyle='dashed')
#	dfe = dfx[~dfx.name.str.startswith('model')].drop('p',axis=1)
#	dfe = pd.pivot_table(dfe,values='y',index='x',aggfunc=np.median)
#	plt.plot(dfe.index,dfe.y,color='black',alpha=0.8,linestyle='dashed')
#	plt.savefig(f'{endval}.png')
#	plt.gcf().clear()
