module DraftsGeostatistics

using CSV
using DataFrames
using Distances
using Distributions: Normal, pdf, cdf
using GeoStats
using LinearAlgebra
using LocalAnisotropies
using NearestNeighbors
using Optim
using Printf
using Setfield
using StableRNGs
using StaticArrays
using StatsBase
using Transducers
using Unitful
import CoordRefSystems: lentype
import NearestNeighbors: MinkowskiMetric
import LocalAnisotropies: LocalKrigingModel, LocalIDWModel

SpatialData = Union{AbstractGeoTable,AbstractDataFrame}
coord_values(pt, axis) = ustrip(to(centroid(pt))[axis])

include("dh/definitions.jl")
include("dh/compositing.jl")
include("dh/desurvey.jl")
include("dh/mergetables.jl")
include("dh/validations.jl")
include("parse_utils.jl")
include("vario.jl")
include("search.jl")
include("blocking.jl")
include("cverror.jl")
include("model.jl")
include("estimates.jl")
include("nlinear.jl")
include("simul.jl")


export composite,
    AdvBallSearch,
    Collar,    
    DrillHole,
    LeaveHoleOut,
    LocalVariogram,
    Interval,
    Survey,

    back_nscore,
    backflag,
    blocks_iterator,
    categorical_compositing,
    cell_declus_tests,
    dgm,
    drillhole,
    exportwarns,
    extract_contacts,
    extract_intrusion_pts,
    extrapolate_borders,
    get_filters,
    global_mean_validation,
    impute_missing_contacts,
    intrusion_model,
    join_anisotropic_variogram,
    local_cverror,
    localaniso_from_pts,
    merge_subblocks,
    multifit,
    nscore,
    read_dh,
    read_expvario,
    read_vartable,
    regblks_estimation,
    regblks_simulation,
    regblks_to_subblks,
    write_dh,
    write_expvario
end
