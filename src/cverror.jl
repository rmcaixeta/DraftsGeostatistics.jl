
struct LeaveHoleOut{B,L} <: GeoStatsValidation.ErrorMethod
    holeids::B
    loss::L
end

struct GroupFolding <: GeoStatsBase.FoldingMethod
    holeids::AbstractVector
end
  
LeaveHoleOut(holeids; loss=Dict()) = LeaveHoleOut(holeids, GeoStatsValidation.assymbol(loss))

function GeoStatsBase.folds(domain::Domain, method::GroupFolding)
    # partition domain
    holeids = method.holeids
    pred(i, j) = isequal(holeids[i], holeids[j])
    partitioner = IndexPredicatePartition(pred)
    p = partition(domain, partitioner)
    s = indices(p)
    n = length(p)

    function pair(i)
        # source and target subsets
        source = [1:(i - 1); (i + 1):n]
        target = [i]

        # indices within subsets
        sinds = reduce(vcat, s[source])
        tinds = reduce(vcat, s[target])

        sinds, tinds
    end

    (pair(i) for i in 1:n)
end
  
function GeoStatsValidation.cverror(setup::GeoStatsValidation.ErrorSetup, geotable::AbstractGeoTable, method::LeaveHoleOut)
    # uniform weights
    weighting = UniformWeighting()

    # ball folds
    folding = GroupFolding(method.holeids)

    wcv = WeightedValidation(weighting, folding, lambda=1, loss=method.loss)

    cverror(setup, geotable, wcv)
end

const INTERPNEIGHBORS = (:minneighbors, :maxneighbors, :neighborhood, :distance)

function local_cverror(model::M, geotable::AbstractGeoTable, method::GeoStatsValidation.ErrorMethod; kwargs...) where {M<:GeoStatsModels.GeoStatsModel}
    I = any(∈(INTERPNEIGHBORS), keys(kwargs)) ? InterpolateNeighbors : Interpolate
    I = model isa LocalKrigingModel ? LocalInterpolate : I
    setup = GeoStatsValidation.InterpSetup{I,M,typeof(kwargs)}(model, kwargs)
    cverror(setup, geotable, method)
  end

function GeoStatsValidation._prediction(s::GeoStatsValidation.InterpSetup{I}, geotable, f) where {I<:LocalInterpolate}
    sdat = view(geotable, f[1])
    sdom = view(domain(geotable), f[2])
    smod = deepcopy(s.model)
    smod = @set smod.localaniso = view(s.model.localaniso, f[2])
    sdat |> I(sdom; model=smod, s.kwargs...)
end


