
function rounded_str(x)
  0.01 <= x < 0.1 ? @sprintf("%.3f", x) :
  0.1 <= x < 10 ? @sprintf("%.2f", x) :
  10 <= x < 100 ? @sprintf("%.1f", x) : 100 <= x < 1000 ? @sprintf("%.0f", x) : @sprintf("%.2e", x)
end

evalstr(x) = eval(Meta.parse(x))

get_filters(tab::AbstractDataFrame, col::Symbol) = [evalstr("row -> " * ln[col]) for ln in eachrow(tab)]
get_filters(rules::AbstractVector{String}) = Dict(d => evalstr("row -> " * d) for d in rules)

function read_vartable(path::String)
  vartable = CSV.read(path, DataFrame)
  for col in names(vartable)
    if eltype(vartable[!, col]) <: AbstractString
      vartable[!, col] .= replace.(vartable[!, col], "“" => "\"", "”" => "\"")
    end
  end
  sort_cols = [x for x in ("var", "model_domain", "pass") if x in names(vartable)]
  vartable |> Sort(sort_cols...)
end

function coord_table(geotable)
  dom = domain(geotable)
  cnames = [:x, :y, :z]
  points = [centroid(dom, i) for i in 1:nelements(dom)]
  ccolumns = map(1:3) do d
    [ustrip(to(p)[d]) for p in points]
  end
  (; zip(cnames, ccolumns)...)
end



## Parquet
# using DuckDB
# using Parquet2: writefile

# function read_parquet(file_path; cols="*")
# 	con = DBInterface.connect(DuckDB.DB, ":memory:")
# 	tab = DBInterface.execute(con,
# 	"""
# 	SELECT $cols
# 	FROM '$(file_path)'
# 	""") |> skipmissing |> DataFrame
# 	DBInterface.close!(con)
# 	tab
# end

# write_parquet(outfile, tab) = writefile(outfile, tab)
## writefile("outfile.parquet", tab)