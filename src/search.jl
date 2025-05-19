



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
            maxpercategory =
                (maxpercategory..., geotable = geotable[:, keys(maxpercategory)])
        end
    end
    AdvBallSearch(domain(geotable), ball; k, maxpercategory, usesectors, rank_metric)
end

GeoStats.KBallSearch(domain, maxneighbors, neighborhood::AdvBallSearch) = neighborhood

AdvBallSearch(geoms; ball, kwargs...) = AdvBallSearch(geoms, ball; kwargs...)

GeoStats.maxneighbors(method::AdvBallSearch) = method.k

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
    k = minimum([length(inds), length(neighbors), method.k])
    rank_metric = method.rank == :same ? metric(method.ball) : Euclidean()
    dists = [
        evaluate(rank_metric, ustrip.(to(pₒ)), ustrip.(to(centroid(domain, ind)))) for
        ind in inds
    ]
    sorted = sortperm(dists)
    inds = inds[sorted]

    # initialize category and sectors constraints if necessary. `Dict` containers
    # are here created with the information necessary for further filtering.
    categs = initcategories(method, inds)
    sectors = initsectors(method)

    # loop each neighbor candidate
    nneigh = 0
    for (i, d) in zip(inds, dists)
        # check category. if exceeded the category, the neighbor is ignored;
        # otherwise categs[:count][col][cat[col]] is incremented in the end after
        # k and maxpersector checks
        if categs[:use]
            cat, pass = Dict(), false
            tab = method.maxpercategory.geotable
            for col in keys(categs[:max])
                col == :geotable && continue
                cat[col] = tab[i, col]
                categs[:count][col][cat[col]] >= categs[:max][col] && (pass = true)
            end
            pass && continue
        end

        # check sectors; if there are no space for neighbors in the indsector, it's
        # ignored; otherwise sectors[:count][indsector] is incremented in the end
        # after k checks
        if sectors[:use]
            centered = sectors[:rotmat]' * (centroid(domain, i) - pₒ) ./ method.ball.radii
            indsector = getsector(centered, sectors[:n])
            sectors[:count][indsector] >= sectors[:max] && continue
        end

        # add neighbor
        nneigh += 1
        neighbors[nneigh] = i
        distances[nneigh] = d * u

        nneigh == k && break

        # add counters
        sectors[:use] && (sectors[:count][indsector] += 1)
        if categs[:use]
            for col in keys(categs[:max])
                col == :geotable && continue
                categs[:count][col][cat[col]] += 1
            end
        end
    end

    nneigh
end



# initialize categories constraints if necessary. e.g. if `maxpercategory =
# (holeid = 2,)`, dict[:count][:holeid] will have all the possible values of
# holeid, initiliazed to zero for further counting.
function initcategories(meth, inds)
    catgs = meth.maxpercategory
    isnothing(catgs) && return Dict(:use => false)
    tab = view(catgs.geotable, inds)
    catvals = Dict(k => unique(tab[:, k]) for k in keys(catgs) if k != :geotable)
    counter = Dict(k => Dict(zip(v, zeros(Int, length(v)))) for (k, v) in catvals)
    Dict(:use => true, :max => catgs, :count => counter)
end

# initialize sectors constraints if necessary. this will create a counter for
# each sector and a rotation matrix to reverse ellipsoid rotation (if an
# ellipsoid search is used). a rotation noise is also added to simplify filtering.
# it is used here because if data is regular and fall onto axes (e.g. cartesian
# grid data), the sector filtering will return misleading results without it.
function initsectors(meth)
    pars = meth.usesectors
    isnothing(pars) && return Dict(:use => false)
    maxpersector = pars.max
    N = embeddim(meth.domain)
    ns = pars.n
    N == 3 && pars.split && (ns *= 2)
    rotmat = isnothing(meth.ball.rotation) ? I : meth.ball.rotation
    Dict(
        :use => true,
        :count => zeros(Int, ns),
        :rotmat => rotmat,
        :max => maxpersector,
        :n => (pars.n, pars.split),
    )
end

# this function returns at which quadrant/octant the given centered coordinate
# is. instead of making multiple if statements, a container is used to get the
# combination of positive/negative coordinates and assign a code to it.
function getsector(v, n)
    angle = atan(v[2], v[1])

    while angle < 0
        angle += 2π
    end
    while angle >= 2π
        angle -= 2π
    end

    sector = 1 + (Int(floor(angle / (2π / n[1]))) % n[1])
    length(v) == 3 && n[2] && v[3] < 0 && (sector += n[1])
    sector
end
