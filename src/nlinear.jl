const QRANGE = LinRange(0.002,0.998,200)


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

function back_nscore(sims::Ensemble, evar, refd)
    reals = mapreduce(vcat, 1:length(sims)) do i
        vals = sims.reals[i][evar]
        (; evar => back_nscore(vals, refd))
    end
    Ensemble(sims.domain, reals, sims.fetch)
end

function back_nscore(gtab::AbstractGeoTable, evar, refd; as_geotable = true)
    vals = getproperty(gtab, evar)
    as_geotable ? georef((; evar => back_nscore(vals, refd)), domain(gtab)) :
    DataFrame((; evar => back_nscore(vals, refd)))
end

back_nscore(vals::AbstractVector, refd) = TableTransforms.qtransform(vals, Normal(), refd)

function find_valid(vec_of_vecs)
	nan_inf_positions = Int[]

	for (i, vec) in enumerate(vec_of_vecs)
		for (j, value) in enumerate(vec)
			if isnan(value) || isinf(value)
				push!(nan_inf_positions, i)
			end
		end
	end
	nan_inf_positions = unique(nan_inf_positions)
	setdiff(1:length(vec_of_vecs), nan_inf_positions)
end

function avg_normal(vn, n::Normal, cache; pipe=Quantile([:Au_ppm]), qs = QRANGE)
  if var(n) == 0
	getproperty(revert(pipe, georef((; vn => [mean(n)])), cache),vn)
  else
	dummy = georef((; vn => quantile(n,qs)))
	vals = getproperty(revert(pipe, dummy, cache),vn)
	mean(vals)
  end
end
avg_normal(vn, n::Number, cache; pipe=0, qs=0) = 0


function merge_normals(normals, weights)
	means = mean.(normals)
	vars = var.(normals)
	sw = sum(weights)

	if sw == 0
	  Normal(0, 0)
	else
	  m = sum(weights .* means) / sw
	  v = sum(weights .* (vars .+ (means .- m).^2)) / sw

	  #m = mean(means)
	  #v = mean(vars) + mean(means .^2) - m^2
	  
	  Normal(m, v ^ 0.5)
	  end
end


function hermite_polynomial_values(k, datavals)
    h = mapreduce(hcat, datavals) do y
      vals = [1.0, -y]
      for i in 2:k-1
        v = -1/sqrt(i) * y * vals[end] - sqrt((i-1)/i) * vals[end-1]
        push!(vals,v)
      end
      vals
    end
end

function hermite_constants(std_vals, ns_vals, poly)
    nhermites = size(poly, 1)
    N = Normal(0, 1)

    coeffs = mapreduce(vcat, 1:nhermites) do i
        if i == 1
            mean(std_vals)
        else
            h = mapreduce(vcat, 2:length(std_vals)) do j
                (std_vals[j-1] - std_vals[j]) / sqrt(i-1) * poly[i-1, j] * pdf(N, ns_vals[j])
            end
            sum(h)
        end
    end
    coeffs
end


function find_R(data_var, f, coeffs; abs_tol=0.00001)
    variance = data_var * f
    sq_coeffs = [i^2 for i in coeffs]
    exponent = [2*i for i in 1:length(coeffs)]
    
    objective(r, exponent, sq_coeffs, variance) = begin
        vector = [r^exponent[i] for i in 1:length(exponent)]
        dot_product = dot(vector, sq_coeffs)
        abs(dot_product - variance)
    end

    opt = optimize(r -> objective(r, exponent, sq_coeffs, variance), -1.0, 1.0, abs_tol=abs_tol)
    minimizer(opt)[1]
end

function transform_dist(r, poly, coeffs)
    changed_constants = [r ^ (2*(p-1)) * coeffs[p] for p in 1:length(coeffs)]
    sum(changed_constants .* poly, dims=1)[:]
end

function dgm(vals, ns_vals, f; nhermites=100, abs_tol=0.00001)
  std_vals = (vals .- mean(vals)) ./ std(vals)
  poly = hermite_polynomial_values(nhermites, ns_vals)
  coeffs = hermite_constants(std_vals, ns_vals, poly)
  r = find_R(var(std_vals), f, coeffs; abs_tol)
  out = transform_dist(r, poly, coeffs)
  std(vals) .* out .+ mean(vals)
end

function dgm(vals, f; kwargs...)
  ns_vals = georef((data=vals,)) |> Quantile(:data)
  ns_vals = ns_vals.data
  dgm(vals, ns_vals, f; kwargs...)
end

function quantiles_dgm(vn, n::Normal, cache, f; pipe=Quantile([:Au_ppm]), qs = QRANGE)
  if var(n) == 0
    getproperty(revert(pipe, georef((; vn => [mean(n) for i in qs])), cache),vn)
  else
    dummy = georef((; vn => quantile(n,qs)))
    vals = getproperty(revert(pipe, dummy, cache),vn)
    dgm(vals, quantile(Normal(0,1),qs), f; nhermites=100, abs_tol=0.00001)
  end
end

function avg_normal_dgm(vn, n::Normal, cache, f; pipe=Quantile([:Au_ppm]), qs = QRANGE)
  dvals = quantiles_dgm(vn, n, cache, f; pipe, qs)
  mean(dvals)
end
avg_normal_dgm(vn, n::Number, cache, f; pipe=0, qs=0) = 0