const QRANGE = LinRange(0.002, 0.998, 200)

struct GMM_pars
  w::AbstractArray
  μ::AbstractArray
  Σ::AbstractArray
end

# wgts = weight(comps, BlockWeighting(55,55,14))
function nscore(gtab::AbstractGeoTable, evar; weights=nothing)
  data = getproperty(gtab, evar)
  wgts = weights isa Union{String,Symbol} ? getproperty(gtab, Symbol(weights)) : weights
  wgts = isnothing(wgts) || wgts isa GeoWeights ? wgts : GeoWeights(domain(gtab), wgts)
  outvec, refd = nscore(data; weights=wgts)
  outtab = georef((; evar => outvec), gtab.geometry)
  outtab, refd
end

function nscore(data; weights=nothing)
  so = TableTransforms.qsmooth(data)
  ss = isnothing(weights) ? so : TableTransforms.qsmooth(quantile(data, weights, LinRange(0, 1, length(data))))
  refd = TableTransforms.EmpiricalDistribution(ss)
  modqtransform(so, refd, Normal()), refd
end

nscore(refd::TableTransforms.EmpiricalDistribution, data) = modqtransform(data, refd, Normal())

function back_nscore(sims::Ensemble, evar, refd)
  reals = mapreduce(vcat, 1:length(sims)) do i
    vals = sims.reals[i][evar]
    (; evar => back_nscore(vals, refd))
  end
  Ensemble(sims.domain, reals, sims.fetch)
end

function back_nscore(gtab::AbstractGeoTable, evar, refd; as_geotable=true)
  vals = getproperty(gtab, evar)
  as_geotable ? georef((; evar => back_nscore(vals, refd)), domain(gtab)) :
  DataFrame((; evar => back_nscore(vals, refd)))
end

back_nscore(vals::AbstractVector, refd) = modqtransform(vals, Normal(), refd)

function modqtransform(values, origin, target)
  # avoid evaluating the quantile at 0 or 1
  nv = length(values)
  smin, smax = LinRange(0.0, 1.0, nv + 2)[[2, nv + 1]]
  pmin = min(0.0 + 1e-3, smin)
  pmax = max(1.0 - 1e-3, smax)
  map(values) do sample
    prob = cdf(origin, sample)
    quantile(target, clamp(prob, pmin, pmax))
  end
end

function avg_normal(n::Normal, refd; qs=QRANGE)
  input = var(n) == 0 ? [mean(n)] : quantile(n, qs)
  mean(back_nscore(input, refd))
end

function avg_normal(n::Normal, dt, gmm::Union{GMM_pars,AbstractGeoTable}, refd; qs=QRANGE)
  bt1 = var(n) == 0 ? dt_backward([dt], [mean(n)], gmm) : dt_backward([dt for i in qs], quantile(n, qs), gmm)
  mean(back_nscore(bt1, refd))
end

function merge_normals(normals, weights)
  means = mean.(normals)
  vars = var.(normals)
  sw = sum(weights)

  if sw == 0
    Normal(0, 0)
  else
    m = sum(weights .* means) / sw
    v = sum(weights .* (vars .+ (means .- m) .^ 2)) / sw

    Normal(m, v^0.5)
  end
end

function hermite_polynomial_values(nk, datavals)
  H = ones(nk + 1, length(datavals))
  H[2, :] .= -1 * datavals

  for k in 3:(nk + 1)
    H[k, :] .= -1 / sqrt(k - 1) .* datavals .* H[k - 1, :] .- sqrt((k - 2) / (k - 1)) .* H[k - 2, :]
  end
  H
end

function hermite_constants(vals, ns_vals, poly)
  nhermites = size(poly, 1)
  N = Normal(0, 1)

  coeffs = mapreduce(vcat, 1:nhermites) do i
    if i == 1
      mean(vals)
    else
      h = mapreduce(vcat, 2:length(vals)) do j
        (vals[j - 1] - vals[j]) / sqrt(i - 1) * poly[i - 1, j] * pdf(N, ns_vals[j])
      end
      sum(h)
    end
  end
  coeffs
end

function find_R(data_var::Float64, f::Float64, coeffs; abs_tol=0.00001)
  variance = data_var * f
  sq_coeffs = [i^2 for i in coeffs[2:end]]
  exponent = [2 * i for i in 1:length(coeffs[2:end])]

  objective(r, exponent, sq_coeffs, variance) = begin
    vector = [r^e for e in exponent]
    dot_product = dot(vector, sq_coeffs)
    abs(dot_product - variance)
  end

  opt = optimize(r -> objective(r, exponent, sq_coeffs, variance), -1.0, 1.0, abs_tol=abs_tol)
  Optim.minimizer(opt)[1]
end

function find_R(vals::AbstractVector, ns_vals::AbstractVector, f; nhermites=100, abs_tol=0.00001, verbose=false)
  poly = hermite_polynomial_values(nhermites, ns_vals)
  coeffs = hermite_constants(vals, ns_vals, poly)
  var_coeffs = sum([x^2 for x in coeffs[2:end]])
  verbose && println("Variance from coeffs: $var_coeffs")

  fvals = f isa AbstractVector ? f : [f]
  mapreduce(vcat, fvals) do f_i
    find_R(var(vals), f_i, coeffs; abs_tol)
  end
end

function transform_dist(r, poly, coeffs)
  changed_constants = [r^(p - 1) * coeffs[p] for p in 1:length(coeffs)]
  sort(sum(changed_constants .* poly, dims=1)[:])
end

function dgm(vals, ns_vals, f; nhermites=100, abs_tol=0.00001)
  poly = hermite_polynomial_values(nhermites, ns_vals)
  coeffs = hermite_constants(vals, ns_vals, poly)
  r = find_R(var(vals), f, coeffs; abs_tol)
  transform_dist(r, poly, coeffs)
end

function dgm(vals, f; kwargs...)
  svals = sort(vals)
  nsvals = LinRange(0.0, 1.0, length(svals) + 2)
  nsvals = quantile(Normal(), nsvals)[2:(length(svals) + 1)]
  dgm(svals, nsvals, f; kwargs...)
end

function quantiles_dgm(n::Normal, refd, f; qs=QRANGE)
  if var(n) == 0
    back_nscore([mean(n) for i in qs], refd)
  else
    vals = back_nscore(quantile(n, qs), refd)
    dgm(vals, f; nhermites=100, abs_tol=0.00001)
  end
end

function quantiles_dgm(n::Normal, dt, gmm::Union{GMM_pars,AbstractGeoTable}, refd, f; qs=QRANGE)
  if var(n) == 0
    bt1 = dt_backward([mean(dt) for i in qs], [mean(n) for i in qs], gmm)
    back_nscore(bt1, refd)
  else
    bt1 = dt_backward(sample(dt, length(qs), replace=true), quantile(n, qs), gmm)
    vals = back_nscore(bt1, refd)
    dgm(vals, f; nhermites=100, abs_tol=0.00001)
  end
end

function fval_la(block_size, γorig, lp, discretization=(4, 4, 4))
  origin, finish = Point(0, 0, 0), Point(block_size...)
  subb = [b / d for (b, d) in zip(block_size, discretization)]
  dgrid = CartesianGrid(origin, finish, Tuple(subb))
  fval_la(dgrid, γorig, lp)
end

function fval_la(dgrid::Domain, γorig, lp)
  samps = centroid.(dgrid)
  γl = [mw_estimator(nothing, γorig, localpair(lp, i)).fun for i in 1:nvals(lp)]
  fvals = mapreduce(vcat, γl) do γi
    avgvar = [γi(p1, p2) for p1 in samps, p2 in samps]
    sill(γorig) - mean(avgvar)
  end
  mean(fvals)
end

function sample_gmm(gmm::GMM_pars, n::Int)
  K, d = size(gmm.μ)
  priors = Distributions.Categorical(gmm.w)
  mvns = [MvNormal(gmm.μ[k, :], Symmetric(gmm.Σ[k])) for k in 1:K]
  mapreduce(x -> rand(mvns[rand(priors)]), hcat, 1:n)
end

function bic(logL, data, ncompo)
  n, d = size(data)
  n_params = (ncompo - 1) + ncompo * d + ncompo * d * (d + 1) ÷ 2
  -2 * logL + n_params * log(n)
end

function gmm_conditional(gmm::GMM_pars, x; first_feature::Bool=true)
  weights, means, covs = (gmm.w, gmm.μ, gmm.Σ)
  n = length(weights)
  cond_means = Vector{Float64}(undef, n)
  cond_stds = Vector{Float64}(undef, n)
  numerators = Vector{Float64}(undef, n)

  for i in 1:length(weights)
    w = weights[i]
    m = means[i, :]
    C = covs[i]

    mu1, mu2 = m[1], m[2]
    sigma11, sigma22, sigma12 = C[1, 1], C[2, 2], C[1, 2]

    if first_feature
      cond_means[i] = mu2 + sigma12 / sigma11 * (x - mu1)
      cond_stds[i] = sqrt(sigma22 - sigma12^2 / sigma11)
      numerators[i] = w * pdf(Normal(mu1, sqrt(sigma11)), x)
    else
      cond_means[i] = mu1 + sigma12 / sigma22 * (x - mu2)
      cond_stds[i] = sqrt(sigma11 - sigma12^2 / sigma22)
      numerators[i] = w * pdf(Normal(mu2, sqrt(sigma22)), x)
    end
  end

  cond_weights = numerators ./ sum(numerators)
  cond_weights, cond_means, cond_stds
end

function conditional_cdf(x2, weights, means, stds)
  mapreduce(+, 1:length(weights)) do i
    weights[i] * cdf(Normal(means[i], stds[i]), x2)
  end
end

function dt_forward(x1_arr, x2_arr, gmm::GMM_pars; threads=1)
  partitions = nreal_partitions(length(x1_arr), threads)
  tmapreduce(vcat, partitions) do inds
    mapreduce(vcat, zip(x1_arr[inds], x2_arr[inds])) do (x1, x2)
      w, m, s = gmm_conditional(gmm, x1)
      cdf2 = conditional_cdf(x2, w, m, s)
      quantile(Normal(), cdf2)
    end
  end
end

function dt_backward(y1_arr, y2_arr, gmm::GMM_pars; bounds=(-5.0, 5.0), method=Roots.Bisection(), threads=1)
  # method Roots.Brent() possible too if preferred
  partitions = nreal_partitions(length(y1_arr), threads)
  tmapreduce(vcat, partitions) do inds
    mapreduce(vcat, zip(y1_arr[inds], y2_arr[inds])) do (y1, y2)
      p = cdf(Normal(), y2)
      w, m, s = gmm_conditional(gmm, y1)
      f(x) = conditional_cdf(x, w, m, s) - p
      if isnothing(bounds)
        μ_mix = sum(w .* m)
        σ_mix = sqrt(sum(w .* (s .^ 2 .+ (m .- μ_mix) .^ 2)))
        fzero(f, quantile(Normal(μ_mix, σ_mix), p))
      else
        find_zero(f, bounds, method)
      end
    end
  end
end

function dt_backward(y1_arr, y2_arr, ref::AbstractGeoTable; power=1.8, neighs=4)
  dom = georef((; :y1 => y1_arr, :y2 => y2_arr), (:y1, :y2))
  pred = ref |> InterpolateNeighbors(domain(dom), model=IDW(power), maxneighbors=neighs) |> values
  pred |> values |> first
end

function dt_table(gmm::GMM_pars; bounds=(-4.5, 4.5), vals=500)
  yi = LinRange(bounds..., vals)
  yij = Iterators.product(yi, yi)
  tcat = (a, b) -> (vcat(a[1], b[1]), vcat(a[2], b[2]))
  y1, y2 = reduce(tcat, yij)
  x2 = dt_backward(y1, y2, gmm, bounds=nothing)
  georef((; :y1 => y1, :y2 => y2, :x2 => x2), (:y1, :y2))
end
