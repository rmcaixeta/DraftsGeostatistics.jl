module DraftsGeostatistics

using CSV
using DataFrames
using Distances
using Distributions: Normal, pdf, cdf
using DuckDB
using GeoStats
using LinearAlgebra
using LocalAnisotropies
using NearestNeighbors
using Optim
using Parquet2: writefile
using Printf
using Setfield
using StableRNGs
using StaticArrays
using StatsBase
using Transducers
using Unitful
import CoordRefSystems: lentype
import NearestNeighbors: MinkowskiMetric
import LocalAnisotropies: LocalKrigingModel, LocalIDWModel, mw_estimator, localpair, nvals

SpatialData = Union{AbstractGeoTable,AbstractDataFrame}
coord_values(pt, axis) = ustrip(to(centroid(pt))[axis])

include("dh/definitions.jl")
include("dh/compositing.jl")
include("dh/desurvey.jl")
include("dh/mergetables.jl")
include("dh/validations.jl")
include("utils.jl")
include("deprecated.jl")
include("vario.jl")
include("search.jl")
include("blocking.jl")
include("cverror.jl")
include("model.jl")
include("estimates.jl")
include("nlinear.jl")
include("simul.jl")
include("declus.jl")

export composite,
  AdvBallSearch,
  Collar,
  DrillHole,
  LeaveHoleOut,
  LocalVariogram,
  Interval,
  Survey,
  avg_normal,
  back_nscore,
  backflag,
  categorical_compositing,
  cell_declus_tests,
  class_blk,
  coord_table,
  dgm,
  drillhole,
  exportwarns,
  extract_contacts,
  extract_intrusion_pts,
  extrapolate_borders,
  fval_la,
  get_filters,
  global_mean_validation,
  impute_missing_contacts,
  intrusion_model,
  join_anisotropic_variogram,
  local_cverror,
  localaniso_from_pts,
  make_subblocks,
  merge_normals,
  merge_subblocks,
  multifit,
  nreal_partitions,
  nscore,
  quantiles_dgm,
  read_dh,
  read_expvario,
  read_table,
  read_vartable,
  #regblks_estimation,
  regblks_simulation,
  regblks_to_subblks,
  variog_ns_to_orig,
  write_dh,
  write_expvario,
  write_parquet
end
