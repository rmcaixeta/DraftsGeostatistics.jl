using DraftsGeostatistics
using GeoStats
using Test

@testset "DraftsGeostatistics.jl" begin

  # DH
  collar_kw = (holeid=:HoleID, x=:Easting, y=:Northing, z=:Elevation, enddepth=:Length_m)
  survey_kw = (holeid=:HoleID, at=:Depth_m, azm=:Azimuth, dip=:Dip)
  litho_kw = (holeid=:holeid, from=:from, to=:to)
  assay_kw = (holeid=:HoleID, from=:From_m, to=:To_m)

  collar = Collar("./data/collar.csv"; collar_kw...)
  survey = Survey("./data/survey.csv"; survey_kw...)
  litho = Interval("./data/litho.csv"; litho_kw...)
  assay = Interval("./data/assay.csv"; assay_kw...)

  dh = drillhole(collar, survey, [litho, assay])

  # MODEL
  dh.table.Code = coalesce.(dh.table.Code, "")
  origin, finish = Point(421965, 7010480, 735), Point(422120, 7010645, 930)
  grid = CartesianGrid(origin, finish, (5, 5, 5))

  int_dict = Dict([("BZPZ", ["BZPZ"]), ("BZUZ", ["BZUZ"])])
  ext_dict = Dict([("BZPZ", ["BZUZ", ""]), ("BZUZ", ["BZPZ", ""])])
  models = Vector{GeoTable}(undef, 2)

  for (i, unit) in enumerate(["BZPZ", "BZUZ"])
    comps = categorical_compositing(
      dh,
      :Code,
      int_dict[unit],
      ext_dict[unit],
      filters=[("Exterior", 1.25), ("Interior", 1.25)]
    )
    sd = extract_intrusion_pts(comps, :info_catg_, ["Interior"], ["Exterior"])
    model = intrusion_model(sd, IDW(1.0), grid, subblocks_split=(2, 2, 2))
    uval = Symbol("$unit")
    model = @transform(model, {uval} = ifelse(ismissing(:SD), -1, Int(:SD < 0)))
    models[i] = (model |> Reject(:SD))
  end

  blks = merge_subblocks(models...)
  blks[!, :Catg] .= "UNK"
  blks[blks.BZUZ .== 1, :Catg] .= "BZUZ"
  blks[blks.BZPZ .== 1, :Catg] .= "BZPZ"
  blks = blks |> Reject(:BZUZ, :BZPZ)

  # EST

end
