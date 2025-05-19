#using GeoStats
#using Distances
#using NearestNeighbors
#import NearestNeighbors: MinkowskiMetric
#import CoordRefSystems: lentype
#using DataFrames


"""
    AdvBallSearch(domain, ball; k=1000, maxpercategory=nothing, usesectors=nothing, rank_metric=:same)
    AdvBallSearch(geotable, ball; k=1000, maxpercategory=nothing, usesectors=nothing, rank_metric=:same)

A method for searching at inside ball (iso- or anisotropic) with multiple possible 
constraints. Constraints available: max `k` neighbors, `maxpercategory` and `usesectors`.
The neighbors are ranked using ball metric if `rank_metric`=:same, otherwise isotropic euclidean.

## Max per category

It can be a `NamedTuple` or `Dict`, where the keys are column name(s) of categorical
properties and the values define the maximum number of neighbors per value of
the property. In the example below, no more than 2 neighbors of same rocktype
are selected. Requires a GeoTable as input.  

`AdvBallSearch(geotable, ball, maxpercategory = (rocktype = 2,))`
```ascii
     _________
  //     . □   \\           ○  Search point
 //      . ▩    \\        ▩ ▲  Neighbors selected
// △ ▲   . ▩     \\       □ △  Neighbors ignored
 ........○........
\\ △ ▲   .       //
 \\      .      //
  \\ △   . □   //
     ‾‾‾‾‾‾‾‾‾
```

## Sectors constraints

Only allow a maximum number of neighbors inside each sector. The sectors 
parameters are: `max` per sector; `n` sectors over ball xy plane; and the 
sectors can be `split` vertically if true - if `n` = 8 and `split` = true, 
a total of 16 sectors will be considered.

`AdvBallSearch(geotable, ball, usesectors=(max=2,n=4,split=false))`
```ascii
     _________
  //     . □   \\           ○  Search point
 //      . ▩    \\        ▩ ▲  Neighbors selected
// ▲ ▲   . ▩     \\       □ △  Neighbors ignored
 ........○........
\\ ▲ ▲   .       //
 \\      .      //
  \\ △   . ▩   //
     ‾‾‾‾‾‾‾‾‾
```
"""
struct AdvBallSearch{D<:Domain,B<:MetricBall,T} <: BoundedNeighborSearchMethod
    # input fields
    domain::D
    k::Int
    ball::B
    maxpercategory::Any
    usesectors::Any
    rank::Symbol
    # state fields
    tree::T
end

function AdvBallSearch(
    domain::D,
    ball::B;
    k::Int = 1000,
    maxpercategory = nothing,
    usesectors = nothing,
    rank_metric = :same,
) where {D<:Domain,B<:MetricBall}
    m = metric(ball)
    xs = [ustrip.(to(centroid(domain, i))) for i = 1:nelements(domain)]
    tree = m isa MinkowskiMetric ? KDTree(xs, m) : BallTree(xs, m)
    AdvBallSearch{D,B,typeof(tree)}(
        domain,
        k,
        ball,
        maxpercategory,
        usesectors,
        rank_metric,
        tree,
    )
end

function AdvBallSearch(
    geotable::GeoTable,
    ball;
    k = 1000,
    maxpercategory = nothing,
    usesectors = nothing,
    rank_metric = :same,
)
    if !isnothing(maxpercategory)
        if !hasproperty(maxpercategory, :geotable)
            catg = keys(maxpercategory)[1]
            maxpercategory = (maxpercategory..., geotable = geotable[:, catg])
        end
    end
    AdvBallSearch(domain(geotable), ball; k, maxpercategory, usesectors, rank_metric)
end

#AdvBallSearch(geoms, ball; kwargs...) = AdvBallSearch(GeometrySet(geoms), ball; kwargs...)

AdvBallSearch(geoms; ball, kwargs...) = AdvBallSearch(geoms, ball; kwargs...)

GeoStats.maxneighbors(method::AdvBallSearch) = method.k

GeoStats.KBallSearch(domain, maxneighbors, neighborhood::AdvBallSearch) = neighborhood


function GeoStats.searchdists!(
    neighbors,
    distances,
    pₒ::Point,
    method::AdvBallSearch;
    mask = nothing,
)
    u = unit(lentype(pₒ))
    domain = method.domain

    inds = search(pₒ, BallSearch(domain, method.ball, method.tree), mask = mask)
    k = minimum([length(inds), length(neighbors)])
    rank_metric = method.rank == :same ? metric(method.ball) : Euclidean()
    cpₒ = ustrip.(to(pₒ))
    dists = [
        evaluate(rank_metric, cpₒ, ustrip.(to(centroid(domain, ind)))) for
        ind in inds
    ]
    sorted = sortperm(dists)
    inds = inds[sorted]
    dists = dists[sorted] .* u

    # apply maxpercategory and maxpersector if necessary
    inds = apply_adv_filters(inds, method, pₒ)

    # loop each neighbor candidate
    nneigh = min(k, length(inds))
    if nneigh > 0
        neighbors[1:nneigh] .= inds[1:nneigh]
        distances[1:nneigh] .= dists[1:nneigh]
    end

    nneigh
end


function apply_adv_filters(inds, meth, pₒ)
    catgs, sectors = meth.maxpercategory, meth.usesectors
    if isnothing(catgs) && isnothing(sectors)
        inds
    elseif isnothing(catgs)
        apply_max_per_sector(inds, meth, pₒ)
    elseif isnothing(sectors)
        apply_max_per_catg(inds, catgs)
    else
        # i1 = apply_max_per_sector(inds, meth, pₒ)
        # i2 = apply_max_per_catg(inds, catgs)
        # intersect(i1, i2)

        i1 = apply_max_per_catg(inds, catgs)
        apply_max_per_sector(i1, meth, pₒ)
    end
end

function apply_max_per_sector(inds, meth, pₒ)
    dom = meth.domain
    radii = meth.ball.radii
    sectors = meth.usesectors
    pars = (sectors.n, sectors.split)
    rotmat = meth.ball.rotation'

    indsector = [rotmat * (centroid(dom, ind) - pₒ) ./ radii for ind in inds]
    indsector = [getsector(ind, pars) for ind in indsector]

    newids = keep_first_n_dups(sectors.max, indsector)
    inds[newids]
end

function apply_max_per_catg(inds, catgs)
    maxval = [catgs[k] for k in keys(catgs) if k != :geotable][1]
    newids = keep_first_n_dups(maxval, view(catgs.geotable, inds))
    inds[newids]
end

function keep_first_n_dups(n, catgs)
    tab = DataFrame((idx = 1:length(catgs), vals = catgs))
    pos = combine(groupby(tab, :vals)) do df
        df[1:min(n, nrow(df)), :]
    end
    sort(pos.idx)
end


# this function returns at which quadrant/octant the given centered coordinate
# is. instead of making multiple if statements, a container is used to get the
# combination of positive/negative coordinates and assign a code to it.
function getsector(v, n)
    angle = atan(v[2], v[1])

    #while angle < 0
    #    angle += 2π
    #end
    #while angle >= 2π
    #    angle -= 2π
    #end
    
    angle = mod(angle, 2π)

    sector = 1 + (Int(floor(angle / (2π / n[1]))) % n[1])
    length(v) == 3 && n[2] && v[3] < 0 && (sector += n[1])
    sector
end
