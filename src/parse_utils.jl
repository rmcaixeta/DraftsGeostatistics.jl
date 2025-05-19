
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
