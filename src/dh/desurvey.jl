# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  	drillhole(collar::Collar, survey::Survey, intervals::Intervals)

Desurvey drill hole based on collar, survey and interval table(s) information.
The intervals can be passed as a single `Interval` or as an array of
`Interval`. Outputs a `DrillHole` object.
"""
function drillhole(collar_::Collar, survey::Survey, intervals::Intervals)
  # pre process information
  collar = deepcopy(collar_)
  pars = getcolnames(survey, intervals)
  warns = validations(collar, survey, intervals)
  in("Error", warns[:, :TYPE]) && (return DrillHole(nothing, nothing, pars, warns))

  # create trace information
  trace = gettrace(collar, survey)
  fillxyz!(trace, pars)

  # merge interval tables
  table = mergetables(intervals, pars)
  fillxyz!(table, trace, pars)

  DrillHole(table, trace, pars, warns)
end

function getcolnames(s, i)
  f = i isa Interval ? i.from : i[1].from
  t = i isa Interval ? i.to : i[1].to
  m = s.method == :tangential
  c = s.convention

  # get most common dip sign and assume it is downwards
  if c == :auto
    df = s.file isa String ? CSV.read(s.file, DataFrame) : s.file
    c = sum(sign.(df[!, s.dip])) > 0 ? :positivedownwards : :negativedownwards
  end

  inv = (c == :positivedownwards)
  pars = (holeid=s.holeid, at=s.at, azm=s.azm, dip=s.dip, from=f, to=t, invdip=inv, tang=m)
end

function gettrace(c, s)
  collar = c.file isa String ? CSV.read(c.file, DataFrame) : c.file
  survey = s.file isa String ? CSV.read(s.file, DataFrame) : s.file

  # rename collar columns to match survey columns if necessary
  n1 = (c.x, c.y, c.z, c.holeid)
  n2 = (:X, :Y, :Z, s.holeid)
  namepairs = [a => b for (a, b) in zip(Symbol.(n1), Symbol.(n2)) if a != b]
  length(namepairs) > 0 && rename!(collar, namepairs...)

  # force coordinates to floats if necessary
  for k in [:X, :Y, :Z]
    force_float = !(eltype(collar[!, k]) <: Float64)
    force_float && (collar[!, k] = convert.(Float64, collar[:, k]))
  end

  # merge collar coordinates to initial survey point
  collar[!, s.at] .= 0.0
  esurv = unique(sort(survey, s.at), s.holeid; keep=:first)
  esurv = esurv[esurv[!, s.at] .> 0, :]
  esurv[!, s.at] .= 0
  survey = vcat(survey, esurv)
  sort!(survey, [s.holeid, s.at])
  trace = leftjoin(survey, collar, on=[s.holeid, s.at])
  sort!(trace, [s.holeid, s.at])
end

# minimum curvature desurvey method
function mincurv(az1, dp1, az2, dp2, d12)
  dp1, dp2 = (90 - dp1), (90 - dp2)

  DL = acos(cosd(dp2 - dp1) - sind(dp1) * sind(dp2) * (1 - cosd(az2 - az1)))
  RF = DL != 0.0 ? 2 * tan(DL / 2) / DL : 1

  dx = 0.5 * d12 * (sind(dp1) * sind(az1) + sind(dp2) * sind(az2)) * RF
  dy = 0.5 * d12 * (sind(dp1) * cosd(az1) + sind(dp2) * cosd(az2)) * RF
  dz = 0.5 * d12 * (cosd(dp1) + cosd(dp2)) * RF
  dx, dy, dz
end

# tangential desurvey method
function tangential(az1, dp1, d12)
  dp1 = (90 - dp1)
  dx = d12 * sind(dp1) * sind(az1)
  dy = d12 * sind(dp1) * cosd(az1)
  dz = d12 * cosd(dp1)
  dx, dy, dz
end

# find survey depths bounding given depth
function findbounds(depths::AbstractArray, at)
  # get closest survey
  nearid = findmin(abs.(depths .- at))[2]
  nearest = depths[nearid]

  # check if depth is after last interval
  nearid == length(depths) && nearest < at && return (nearid, nearid)
  #println("$nearest $at")

  # return (previous, next) survey ids for given depth
  nearest == at && return (nearid, nearid)
  #nearest >  at && return (nearid-1, nearid)
  nearest > at && nearid > 1 && return (nearid - 1, nearid)
  nearest > at && nearid < 1.5 && return (nearid, nearid)
  nearest < at && return (nearid, nearid + 1)
end

# convert survey angles to 3-D vector and vice versa
angs2vec(az, dp) = [sind(az) * cosd(-dp), cosd(az) * cosd(-dp), sind(-dp)]
vec2angs(i, j, k) = [atand(i, j), -asind(k)]

# average angle between two surveyed intervals
function weightedangs(angs1, angs2, d12, d1x)
  # angle to vectors
  v1 = angs2vec(angs1...)
  v2 = angs2vec(angs2...)

  # weight average vector according to distance to surveys
  p2 = d1x / d12
  p1 = 1 - p2
  v12 = v1 * p1 + v2 * p2
  v12 /= sqrt(sum(abs2, v12))

  # convert average vector to survey angles and return it
  azm, dip = vec2angs(v12...)
  azm, dip
end

# fill xyz for dh trace files
function fillxyz!(trace, pars)
  # get column names
  at, az, dp, tang = pars.at, pars.azm, pars.dip, pars.tang
  f = pars.invdip ? -1 : 1

  # loop trace file
  for i in 1:size(trace, 1)
    # pass depth 0 where collar coordinates are already available
    trace[i, at] == 0 && continue

    # get distances and angles; return increments dx,dy,dz
    d12 = trace[i, at] - trace[i - 1, at]
    az1, dp1 = trace[i - 1, az], f * trace[i - 1, dp]
    az2, dp2 = trace[i, az], f * trace[i, dp]
    dx, dy, dz = tang ? tangential(az1, dp1, d12) : mincurv(az1, dp1, az2, dp2, d12)

    # add increments dx,dy,dz to previous coordinates
    trace[i, :X] = dx + trace[i - 1, :X]
    trace[i, :Y] = dy + trace[i - 1, :Y]
    trace[i, :Z] = dz + trace[i - 1, :Z]
  end
end

# fill xyz for interval tables with from-to information
function fillxyz!(tab, trace, pars; output=["mid"])
  # get column names
  bh, at, az, dp, tang = pars.holeid, pars.at, pars.azm, pars.dip, pars.tang
  f = pars.invdip ? -1 : 1
  sufix = Dict("mid" => "", "from" => "_FROM", "to" => "_TO")

  # initialize coordinate columns with float values
  for out in output
    tab[!, Symbol("X" * sufix[out])] .= -9999.9999
    tab[!, Symbol("Y" * sufix[out])] .= -9999.9999
    tab[!, Symbol("Z" * sufix[out])] .= -9999.9999
  end

  # get first hole name and get trace of that hole
  lastbhid = tab[1, bh]
  dht = trace[(trace[!, bh] .== lastbhid), :]

  # loop all intervals
  for i in 1:size(tab, 1)
    # get hole name and mid point depth
    bhid = tab[i, bh]
    # update trace if hole name is different than previous one
    bhid != lastbhid && (dht = trace[(trace[!, bh] .== bhid), :])
    lastbhid = bhid
    # pass if no survey is available (WARN)
    size(dht, 1) == 0 && continue

    # get surveys bounding given depth
    for out in output
      atx =
        out == "mid" ? tab[i, pars.from] + tab[i, :LENGTH] / 2 :
        out == "from" ? tab[i, pars.from] : out == "to" ? tab[i, pars.to] : missing
      ismissing(atx) && error("Missing info; and output kw must be a list containing \"mid\", \"from\" or \"to\"")
      b = findbounds(dht[:, at], atx)
      d1x = atx - dht[b[1], at]

      if d1x == 0
        # if interval depth matches trace depth, get trace coordinates
        tab[i, Symbol("X" * sufix[out])] = dht[b[1], :X]
        tab[i, Symbol("Y" * sufix[out])] = dht[b[1], :Y]
        tab[i, Symbol("Z" * sufix[out])] = dht[b[1], :Z]
      else
        # if not, calculate coordinates increments dx,dy,dz
        d12 = dht[b[2], at] - dht[b[1], at]
        az1, dp1 = dht[b[1], az], f * dht[b[1], dp]
        az2, dp2 = dht[b[2], az], f * dht[b[2], dp]
        azx, dpx = b[1] == b[2] ? (az2, dp2) : weightedangs([az1, dp1], [az2, dp2], d12, d1x)
        dx, dy, dz = tang ? tangential(az1, dp1, d1x) : mincurv(az1, dp1, azx, dpx, d1x)

        # add increments dx,dy,dz to trace coordinates
        tab[i, Symbol("X" * sufix[out])] = dx + dht[b[1], :X]
        tab[i, Symbol("Y" * sufix[out])] = dy + dht[b[1], :Y]
        tab[i, Symbol("Z" * sufix[out])] = dz + dht[b[1], :Z]
      end
    end
  end
  # check if some coordinate was not filled and return a warning if necessary
  filter!(row -> row[Symbol("X" * sufix[output[1]])] != -9999.9999, tab)
end
