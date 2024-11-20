

## MAIN WORKFLOW: categorical_compositing() -> extract_intrusion_pts() -> intrusion_model()

function categorical_compositing(
    dh_,
    var,
    interior,
    exterior;
    filters = [("Exterior", 5.0), ("Interior", 2.1)],
)
    pars = dh_.pars
    bh_, fr_, to_ = pars.holeid, pars.from, pars.to
    dh = dh_.table |> Sort(bh_, fr_)

    # filters
    nvals = nrow(dh)
    interior_filter = findall(x -> x in interior, dh[!, var])
    exterior_filter = findall(x -> x in exterior, dh[!, var])

    catgtype = ["Ignored" for i = 1:nvals]
    catgtype[interior_filter] .= "Interior"
    catgtype[exterior_filter] .= "Exterior"

    ignore_filter = [x for x = 1:nvals if !(x in vcat(interior_filter, exterior_filter))]

    # table grouped by codes
    code = set_interval_code_(dh, catgtype, bh_)
    dht = hcat(dh, DataFrame((code_catg_ = code, info_catg_ = catgtype)))
    grp = groupby(dht, [:code_catg_, :info_catg_, bh_])
    grp = combine(grp, fr_ => minimum => fr_, to_ => maximum => to_)
    grp.LENGTH = grp[!, to_] - grp[!, fr_]

    # return here if ignored is wanted and no filters of intervals

    ## ignore ignored
    grp = grp[grp[!, :info_catg_].!="Ignored", :]

    # get upper and lower bounds of the hole, to not filter
    holeids = unique(grp[!, bh_])
    bounds = vcat(
        [searchsortedfirst(grp[!, bh_], x) for x in holeids],
        [searchsortedlast(grp[!, bh_], x) for x in holeids],
    )

    # apply filters for small intervals
    for (ctype, cleng) in filters

        filter_t = (grp[!, :info_catg_] .== ctype) .& (grp[!, :LENGTH] .< cleng)
        filter_t = [x for x in findall(filter_t) if !(x in bounds)]

        filter_t = [x for x = 1:nrow(grp) if !(x in filter_t)]
        grp = grp[filter_t, :]
        code = set_interval_code_(grp, grp[!, :info_catg_], bh_)
        grp = hcat(select(grp, Not(:code_catg_)), DataFrame((code_catg_ = code,)))

        grp = groupby(grp, [:code_catg_, :info_catg_, bh_])
        grp = combine(grp, fr_ => minimum => fr_, to_ => maximum => to_)
        grp.LENGTH = grp[!, to_] - grp[!, fr_]

    end
    fillxyz!(grp, dh_.trace, pars, output = ["mid", "from", "to"])
    DrillHole(grp, dh_.trace, pars, dh_.warns)
end


function set_interval_code_(dh, catg, bh_)
    # aux function; this loop does not consider gaps currently
    # dh must be already sorted
    code = [1]
    for i = 2:nrow(dh)
        c = code[end]
        if dh[i-1, bh_] == dh[i, bh_] && catg[i-1] == catg[i]
            push!(code, c)
        else
            push!(code, c + 1)
        end
    end
    code
end

function extract_intrusion_pts(
    dh_,
    var,
    interior,
    exterior,
    composite = true,
    clipval = -999.9;
    compleng = 1.0,
)

    dh = copy(dh_.table)
    bh_, fr_ = dh_.pars.holeid, dh_.pars.from

    interior_filter = findall(x -> x in interior, dh[!, var])
    exterior_filter = findall(x -> x in exterior, dh[!, var])
    bhids = dh[!, bh_]

    if composite
        # split intervals in compleng m points or minimum of two values per interval; this rule can be improved
        interior_cmp = mapreduce(vcat, interior_filter) do i
            lngt, bh = abs(dh[i, :LENGTH]), dh[i, bh_]
            fr = dh[i, fr_]
            parts = ceil(Int, lngt / compleng) + 1
            intervals = LinRange(0.01, lngt - 0.01, parts)
            ofrom = [fr + itx for itx in intervals]
            obhid = [bh for itx in intervals]
            tab = DataFrame(bh_ => obhid, fr_ => ofrom)
            fillxyz!(tab, dh_.trace, dh_.pars, output = ["from"])
            [(row[bh_], row[:X_FROM], row[:Y_FROM], row[:Z_FROM]) for row in eachrow(tab)]
        end

        exterior_cmp = mapreduce(vcat, exterior_filter) do i
            lngt, bh = abs(dh[i, :LENGTH]), dh[i, bh_]
            fr = dh[i, fr_]
            parts = ceil(Int, lngt / compleng) + 1
            intervals = LinRange(0.01, lngt - 0.01, parts)
            ofrom = [fr + itx for itx in intervals]
            obhid = [bh for itx in intervals]
            tab = DataFrame(bh_ => obhid, fr_ => ofrom)
            fillxyz!(tab, dh_.trace, dh_.pars, output = ["from"])
            [(row[bh_], row[:X_FROM], row[:Y_FROM], row[:Z_FROM]) for row in eachrow(tab)]
        end

        geoms = vcat(
            [SVector(x[2:end]...) for x in interior_cmp],
            [SVector(x[2:end]...) for x in exterior_cmp],
        )
        bhids = vcat([x[1] for x in interior_cmp], [x[1] for x in exterior_cmp])
        nint, nall = [length(interior_cmp), length(geoms)]
        interior_filter = 1:nint
        exterior_filter = (nint+1):nall
    else
        if "Y_FROM" in names(dh)
            geoms =
                [SVector(row[:X_FROM], row[:Y_FROM], row[:Z_FROM]) for row in eachrow(dh)]
        else
            # if not in columns, run fillxyz first
            fillxyz!(dh, dh_.trace, dh_.pars, output = ["from"])
            geoms =
                [SVector(row[:X_FROM], row[:Y_FROM], row[:Z_FROM]) for row in eachrow(dh)]
        end
    end

    interior = geoms[interior_filter]
    exterior = geoms[exterior_filter]

    # interior signed distances
    tree = KDTree(exterior) #origin
    idxs, dists = nn(tree, interior) #target
    interior_sd = -1 * dists

    # exterior signed distances
    tree = KDTree(interior) #origin
    idxs, dists = nn(tree, exterior) #target
    exterior_sd = dists

    # clip vals
    clipval == -999.9 && (clipval = std(vcat(interior_sd, exterior_sd)) * 2)
    interior_sd = clamp.(interior_sd, -clipval, 0)
    exterior_sd = clamp.(exterior_sd, 0, clipval)

    sd = fill(NaN, length(geoms))
    sd[interior_filter] .= interior_sd
    sd[exterior_filter] .= exterior_sd
    out = DataFrame(
        bh_ => bhids,
        :X => [x[1] for x in geoms],
        :Y => [x[2] for x in geoms],
        :Z => [x[3] for x in geoms],
        :SD => round.(sd, digits = 2),
    )
    georef(out, (:X, :Y, :Z))
end

LocalEstimator = LocalAnisotropies.LocalKrigingModel

function intrusion_model(
    sd,
    interpolant,
    grid;
    subblocks_split = (3, 3, 8),
    sub_interpolant = IDW(3),
    maxneighbors = 100,
    sub_maxneighbors = 50,
)

    sdev = std(sd.SD)
    dh = @transform(sd, :SD = :SD / sdev)
    res_init = get_spacing(grid)

    # do a 1st round estimate in upscaled blocks
    ##

    # estimate in parent blocks
    out = grid isa CartesianGrid ? grid : grid.geometry
    interpfun = interpolant isa LocalEstimator ? LocalInterpolate : InterpolateNeighbors
    est = dh |> interpfun(out, :SD => interpolant, maxneighbors = maxneighbors)
    #add ijk
    est = hcat(est, georef((IJK = 1:nrow(est),), est.geometry))
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
        sub_interpolant_ = LocalKriging(sub_interpolant.method,lp_,sub_interpolant.γ)
        dh |> LocalInterpolate(
            refined.geometry,
            :SD => sub_interpolant_,
            maxneighbors = sub_maxneighbors,
        )
    else
        chunks = collect(
            partition(refined.geometry, UniformPartition(Threads.nthreads(), false)),
        )
        foldxt(
            vcat,
            Transducers.Map(
                x ->
                    dh |> InterpolateNeighbors(
                        x,
                        :SD => sub_interpolant,
                        maxneighbors = sub_maxneighbors,
                    ),
            ),
            chunks,
        )
    end

    est = hcat(est, georef((IJK = refined.IJK,), refined.geometry))
    est = vcat(est, ok)
    est = @transform(est, :SD = :SD * sdev)
    est
end



function extract_contacts(grp; valid_holeids = nothing)#, minz=-8000, maxz=8000)
    ## i guess it works only if the stratigraphy is exterior at the bottom and interior at the top
    pars = grp.pars
    bh_, fr_, to_ = pars.holeid, pars.from, pars.to
    maxz = maximum(vcat(grp.table[!, :Z_FROM], grp.table[!, :Z_TO])) + 1
    minz = minimum(vcat(grp.table[!, :Z_FROM], grp.table[!, :Z_TO])) - 1

    dh_contact = groupby(grp.table, [bh_, :info_catg_])
    dh_contact = combine(
        dh_contact,
        fr_ => minimum => fr_,
        to_ => maximum => to_,
        :Z_FROM => first => :MinZ,
        :Z_TO => last => :MaxZ,
        :Z_FROM => maximum => :Z0,
    )
    dh_contact.MinZ .-= 0.01
    dh_contact.MaxZ .+= 0.01

    ## The folowing might be more correct, but need tests
    #dh_contact = combine(dh_contact, fr_ => minimum => fr_, to_ => maximum => to_, :Z_FROM => minimum => :MinZ_FROM, :Z_FROM => maximum => :MaxZ_FROM, :Z_TO => minimum => :MinZ_TO, :Z_TO => maximum => :MaxZ_TO)
    #dh_contact.Z0 = maximum(df[!, [:MaxZ_FROM, :MaxZ_TO]], dims=2) #[:]
    #dh_contact.MinZ = maximum(df[!, [:MaxZ_FROM, :MaxZ_TO]], dims=2) .- 0.01
    #dh_contact.MaxZ = minimum(df[!, [:MinZ_FROM, :MinZ_TO]], dims=2) .+ 0.01

    interior = dh_contact[!, :info_catg_] .== "Interior"
    dh_contact[interior, :MinZ] .= minz
    dh_contact[.!interior, :MaxZ] .= maxz
    dh_contact = groupby(dh_contact, bh_)
    dh_contact = combine(
        dh_contact,
        fr_ => minimum => fr_,
        to_ => maximum => to_,
        :MinZ => maximum => :MinZ,
        :MaxZ => minimum => :MaxZ,
        :info_catg_ => first => :First,
        :info_catg_ => last => :Last,
        :info_catg_ => length => :COUNT,
        :Z0 => maximum => :Z0,
    )

    contact_rule =
        (dh_contact[!, :COUNT] .> 1) .& (dh_contact[!, :First] .!= dh_contact[!, :Last])
    if valid_holeids != nothing
        valid_holeids = [x in valid_holeids for x in dh_contact[!, bh_]]
        contact_rule = contact_rule .& valid_holeids
    end
    holeids_contact = dh_contact[contact_rule, bh_]
    holeids_nocontact = dh_contact[.!contact_rule, bh_]

    # Get Z contact; add some rule later to get first contact or last contact
    contact_pts = filter(row -> row[bh_] in holeids_contact, grp.table)
    contact_pts = groupby(contact_pts, [bh_, :info_catg_])
    contact_pts = combine(contact_pts, :code_catg_ => last => :code)

    contact_pts = contact_pts[contact_pts[!, :info_catg_].=="Exterior", :code]
    contact_pts = filter(row -> row.code_catg_ in contact_pts, grp.table)
    contact_pts = select(contact_pts, bh_, :X_FROM => :X, :Y_FROM => :Y, :Z_FROM => :Z)

    # Merge info
    dh_contact = leftjoin(dh_contact, contact_pts, on = bh_)
    notmiss = .!ismissing.(dh_contact.Z)
    dh_contact[notmiss, :MinZ] .= dh_contact[notmiss, :Z]
    dh_contact[notmiss, :MaxZ] .= dh_contact[notmiss, :Z]
    badrule = dh_contact.MinZ .> dh_contact.MaxZ
    dh_contact[badrule, :MinZ] .= minz
    dh_contact[badrule, :MaxZ] .= maxz
    select!(dh_contact, Not(:COUNT))
    dh_contact
end


function impute_missing_contacts(contact_table, interpolant, comps; only_2d = false)
    pars = comps.pars
    bh_, fr_, to_, dp_ = pars.holeid, pars.from, pars.to, pars.dip
    maxz = maximum(contact_table.MaxZ)

    no_contact = filter(row -> ismissing(row.Z), contact_table)[!, bh_]
    to_estim = combine(
        groupby(comps.trace, bh_),
        dp_ => mean => dp_,
        :X => mean => :X,
        :Y => mean => :Y,
    )
    to_estim = to_estim[abs.(to_estim[!, dp_]).==90, [bh_, :X, :Y]]
    to_estim = filter(row -> row[bh_] in no_contact, to_estim)
    hd = filter(row -> !ismissing(row.Z), contact_table)

    (nrow(to_estim) == 0 || nrow(hd) == 0) && (return comps)

    to_estim = leftjoin(to_estim, select(contact_table, Not(:X, :Y, :Z)), on = bh_)
    to_estim = georef(to_estim, (:X, :Y))
    hd = georef(hd, (:X, :Y))

    est =
        hd |>
        UniqueCoords(:Z => mean) |>
        InterpolateNeighbors(to_estim.geometry, :Z => interpolant, maxneighbors = 100)
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

    newdh = vcat(comps.table, newrows, cols = :intersect) |> Sort([bh_, fr_])
    newdh[!, :code_catg_] = set_interval_code_(newdh, newdh[!, :info_catg_], bh_)
    newdh = groupby(newdh, [:code_catg_, :info_catg_, bh_])
    newdh = combine(newdh, fr_ => minimum => fr_, to_ => maximum => to_)
    newdh[!, :LENGTH] = newdh[!, to_] - newdh[!, fr_]

    fillxyz!(newdh, comps.trace, pars, output = ["mid", "from", "to"])
    DrillHole(newdh, comps.trace, pars, comps.warns)

end

# length to extrapolation assumes it is subvertical, does not use survey angle to adjust it
function extrapolate_borders(dh_; min_dip = 45.0, ignore_lower = [], ignore_upper = [])
    pars = dh_.pars
    bh_, fr_, to_, dp_ = pars.holeid, pars.from, pars.to, pars.dip
    dh = dh_.table |> Sort([bh_, fr_])
    maxz = maximum(dh.Z_FROM)
    minz = minimum(dh.Z_TO)

    dh_valid = combine(groupby(dh_.trace, bh_), dp_ => first => dp_)
    dh_valid = dh_valid[abs.(dh_valid[!, dp_]).>min_dip, bh_]
    dh_valid = filter(row -> row[bh_] in dh_valid && !(row[bh_] in ignore_upper), dh)

    upper_extrap = groupby(dh_valid, bh_)
    upper_extrap = combine(
        upper_extrap,
        :code_catg_ => first => :code_catg_,
        :info_catg_ => first => :info_catg_,
        fr_ => first => to_,
        :Z_FROM => first => :Z,
    )
    upper_extrap[!, :LENGTH] = maxz .- upper_extrap[!, :Z]
    upper_extrap[!, fr_] = upper_extrap[!, to_] - upper_extrap[!, :LENGTH]
    upper_extrap = upper_extrap[upper_extrap[!, :info_catg_].=="Interior", Not(:Z)]

    dh_valid = combine(groupby(dh_.trace, bh_), dp_ => last => dp_)
    dh_valid = dh_valid[abs.(dh_valid[!, dp_]).>min_dip, bh_]
    dh_valid = filter(row -> row[bh_] in dh_valid && !(row[bh_] in ignore_lower), dh)

    lower_extrap = groupby(dh_valid, bh_)
    lower_extrap = combine(
        lower_extrap,
        :code_catg_ => last => :code_catg_,
        :info_catg_ => last => :info_catg_,
        to_ => last => fr_,
        :Z_TO => last => :Z,
    )
    lower_extrap[!, :LENGTH] = lower_extrap[!, :Z] .- minz
    lower_extrap[!, to_] = lower_extrap[!, fr_] + lower_extrap[!, :LENGTH]
    lower_extrap = lower_extrap[lower_extrap[!, :info_catg_].=="Exterior", Not(:Z)]

    newdh = vcat(dh, upper_extrap, lower_extrap, cols = :intersect) |> Sort([bh_, fr_])
    newdh[!, :code_catg_] = set_interval_code_(newdh, newdh[!, :info_catg_], bh_)
    newdh = groupby(newdh, [:code_catg_, :info_catg_, bh_])
    newdh = combine(newdh, fr_ => minimum => fr_, to_ => maximum => to_)
    newdh[!, :LENGTH] = newdh[!, to_] - newdh[!, fr_]

    fillxyz!(newdh, dh_.trace, pars, output = ["mid", "from", "to"])
    DrillHole(newdh, dh_.trace, pars, dh_.warns)
end
