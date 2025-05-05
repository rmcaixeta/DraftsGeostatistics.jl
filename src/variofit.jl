

"""
    multifit(structures, experimental_vario; sill=1.0)

An extension of `fit` method that accept multiple structures (aka `NestedVariograms`)
or multiple variograms of a same variable that should share same nugget effect and
sill (ie. variogram of Au at two different directions). `structures` and
`experimental_vario` should be informed as array

Example:

e1 = DirectionalVariogram((1.,0.), dh, :Au, maxlag=20,nlags=10)
e2 = DirectionalVariogram((0.,1.), dh, :Au, maxlag=20,nlags=10)
structures = [SphericalVariogram, ExponentialVariogram]
fitted_variograms = multifit(structures, [e1,e2])
"""
function multifit(sts::AbstractVector, evario::AbstractVector; sill::Number = 1.0, fixnugget=-1.0)
    initguess = initialguess(sts, evario) #[cc1, cc2, r1_e1, r2_e1, r1_e2, r2_e2]
    opt = optimize(x -> obj_function(x, evario, sts, sill, fixnugget), initguess)
    res = Optim.minimizer(opt)
    [model_from_pars(res, sill, sts; id = x, fixnugget)[1] for x = 1:length(evario)]
end

function model_from_pars(pars, vvar, sts; id = 1, fixnugget=-1.0)
    nst = length(sts)
    cc = pars[1:nst]
    rang = pars[nst*id+1:nst*(id+1)]
    ngt = fixnugget < 0 ? 1 - sum(cc) : fixnugget
    cc = fixnugget < 0 ? cc : (cc ./ sum(cc)) .* (1-ngt)
    model = (ngt * vvar) * NuggetEffect()
    for i = 1:nst
        model += sts[i](sill = cc[i] * vvar, range = rang[i])
    end
    flagerror = (sum(pars .< 0) > 0) | (sum(cc) > 1.0)
    model, flagerror
end

variovalues(v::EmpiricalVariogram) = (v.abscissas, v.ordinates, v.counts)

function obj_function(pars, vario, sts, vvar, fixnugget)
    mse = 0.0
    for (i, v) in enumerate(vario)
        model, flag = model_from_pars(pars, vvar, sts; id = i, fixnugget)
        penalty = flag ? vvar*100 : 0.0

        x, y, n = variovalues(v)
        x = x[n.>0]
        y = y[n.>0]
        n = n[n.>0]

        error = (model.(x) .- y) .^ 2
        wgt = n .* (1.0 ./ x)
        mse += sum(error .* wgt) / sum(wgt)
        mse += penalty
    end
    mse
end

function initialguess(sts, evario)
    modl = [GeoStatsFunctions.fit(sts[1], ev, h -> exp(-ustrip(h) / 100)) for ev in evario]
    pars = [0.9 / length(sts) for i in sts] #cc
    for m in modl
        for s = 1:length(sts)
            push!(pars, s * ustrip(range(m)) / length(sts))
        end
    end
    pars
end


function join_into_aniso_variogram(γ, anisotropy_pattern)

    ngt = "$(rounded_str(nugget(γ[1])))*NuggetEffect()"
    p = structures(γ[1])
    cc = [rounded_str(x * y.sill) for (x, y) in zip(p[2], p[3])]
    sts = [string(GeoStatsFunctions.constructor(x)) for x in p[3]]
    sts = [replace(x, "GeoStatsFunctions." => "") for x in sts]
    ball = [anisotropy_pattern for x in p[3]]

    for (i, v) in enumerate(γ)
        ranges = Tuple(ustrip(range(x)) for x in structures(v)[3])
        for (j, r) in enumerate(ranges)
            ball[j] = replace(ball[j], "evario$i" => rounded_str(r))
        end
    end

    #out = ngt + sum(c*v(evalstr(b)) for (c, v, b) in zip(cc, sts, ball))
    sts = ["$c*$v($b)" for (c, v, b) in zip(cc, sts, ball)]
    insert!(sts, 1, ngt)
    join(sts, " + ")
end


struct VariographyResult
    groups::AbstractVector{AbstractString}
    exp::AbstractVector{AbstractVector{EmpiricalVariogram}}
    fit::AbstractVector{AbstractVector{Variogram}}
    out::AbstractDataFrame
end

Base.length(x::VariographyResult) = length(x.groups)

function multifit(dh::AbstractGeoTable, vartable::AbstractDataFrame; ns = false)
    filters = get_filters(vartable, :domain)
    multifit(dh, vartable, filters; ns)
end

function multifit(
    dh::AbstractGeoTable,
    vartable::AbstractDataFrame,
    filters::AbstractVector;
    ns = false,
    localpars = nothing,
    localpts = nothing,
)
    enames = [v for v in names(vartable) if startswith(v, "evario")] #sort or constrain?

    vname = Vector{String}()
    vfit = Vector{Vector{Variogram}}()
    vexp = Vector{Vector{EmpiricalVariogram}}()
    vfun = Vector{String}()
    vtable = vartable[:, [:var, :domain]]

    for (ix, ln) in enumerate(eachrow(vartable))
        #fx = evalstr("row->" * ln.domain)
        dhx = dh |> GeoStats.Filter(filters[ix]) |> Select(ln.var) |> DropMissing()
        ns && (dhx = dhx |> Quantile())
        sts = evalstr(ln.fitting_structs)
        evario = mapreduce(vcat, enames) do col
            pars = evalstr(ln[col])
            args1 = haskey(pars, :args) ? pars.args : []
            args2 = (dhx, ln.var) # modify for multivar case
            kws2 = pars.fun == LocalVariogram ? (localpars=localpars, localpts=localpts) : []
            pars.fun(args1..., args2...; kws2..., pars.kwargs...)
        end
        !(evario isa AbstractVector) && (evario = [evario])
        γ = multifit(sts, evario; sill = var(getproperty(dhx, ln.var)))
        push!(vname, ln.var * "_" * ln.domain)
        push!(vfit, γ)
        push!(vexp, evario)
        push!(vfun, join_into_aniso_variogram(γ, ln.anisotropy))
    end

    vtable[!, :variomodel] = vfun
    VariographyResult(vname, vexp, vfit, vtable)
end

function LocalVariogram(dhx, varn; localpars, localpts, kwargs...)
    spars = nnpars(localpars, localpts, dhx)
    localvariography(dhx, spars, varn; kwargs...)
end

function join_anisotropic_variogram(γs, rotation_3d)
    ngt, cc, mods = structures(γs[1])
    mods = [typeof(m).name.wrapper for m in mods]
    rang = mapreduce(vcat, γs) do γ
        p = structures(γ)
        mapreduce(γi -> range(γi), hcat, p[3])
    end
    NuggetEffect(ngt) + sum(cc[i] * mods[i](ranges = Tuple(vec(rang[:,i])), rotation = rotation_3d) for i in 1:length(cc))
end