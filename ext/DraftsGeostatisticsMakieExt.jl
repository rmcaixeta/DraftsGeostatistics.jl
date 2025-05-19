# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DraftsGeostatisticsMakieExt

using DraftsGeostatistics
import Makie

SD = DraftsGeostatistics.SpatialData
cval = DraftsGeostatistics.coord_values

Makie.@recipe(PLT, obj, section) do scene
  return Makie.Attributes(; alpha=0.6, size=5.0, color=:red, colormap=:jet, buffer=5, space=:data)
end

Makie.plottype(::SD, ::String) = PLT{<:Tuple{SD,String}}

function get_filter(obj, ax, lims)
  vals = if hasproperty(obj, Symbol(ax))
    getproperty(obj, Symbol(ax))
  elseif hasproperty(obj, Symbol(lowercase(ax)))
    getproperty(obj, Symbol(lowercase(ax)))
  else
    axid = Dict("X" => 1, "Y" => 2, "Z" => 3)[ax]
    [cval(x, axid) for x in obj.geometry]
  end
  lims[1] .<= vals .<= lims[2]
end

function get_proj(obj, ax, f)
  u, v = Dict("X" => "YZ", "Y" => "XZ", "Z" => "XY")[ax]
  if hasproperty(obj, Symbol(ax))
    getproperty(obj, Symbol(u))[f], getproperty(obj, Symbol(v))[f]
  elseif hasproperty(obj, Symbol(lowercase(ax)))
    ul, vl = lowercase.((u, v))
    getproperty(obj, Symbol(ul))[f], getproperty(obj, Symbol(vl))[f]
  else
    iaxis = Dict('X' => 1, 'Y' => 2, 'Z' => 3)
    iu, iv = iaxis[u], iaxis[v]
    pts = obj[f, :geometry]
    [cval(x, iu) for x in pts], [cval(x, iv) for x in pts]
  end
end

function get_color(obj, c, f)
  if hasproperty(obj, Symbol(c))
    Float64.(getproperty(obj, Symbol(c))[f])
  else
    c
  end
end

function preprocess(obj, section, buffer, color)
  sect = split(replace(section, r"\s+" => ""), "=")
  axis = uppercase(sect[1])
  sect = parse(Float64, sect[end])
  lims_md = (sect - buffer/2, sect + buffer/2)
  f_md = get_filter(obj, axis, lims_md)
  xm, ym = get_proj(obj, axis, f_md)
  c_md = get_color(obj, color, f_md)

  xm, ym, c_md
end

function Makie.convert_arguments(P::Type{<:PLT}, obj::SD, section::String)
  obj, section
end

function Makie.plot!(plot::PLT)
  # input pars
  obj = plot[1]
  section = plot[2]
  buffer = plot[:buffer]
  color = plot[:color]

  out = Makie.@lift preprocess($obj, $section, $buffer, $color)
  x = Makie.@lift $out[1]
  y = Makie.@lift $out[2]
  c = Makie.@lift $out[3]

  Makie.scatter!(
    plot,
    x,
    y,
    color=c,
    markersize=plot[:size],
    markerspace=plot[:space],
    alpha=plot[:alpha],
    colormap=plot[:colormap]
  )
end

end
