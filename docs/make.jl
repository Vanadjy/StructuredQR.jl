using StructuredQR
using Documenter

DocMeta.setdocmeta!(StructuredQR, :DocTestSetup, :(using StructuredQR); recursive = true)

makedocs(;
  modules = [StructuredQR],
  doctest = true,
  linkcheck = false,
  strict = false,
  authors = "Valentin Dijon <valentin.dijon@polymtl.ca>",
  repo = "https://github.com/JuliaSmoothOptimizers/StructuredQR.jl/blob/{commit}{path}#{line}",
  sitename = "StructuredQR.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://JuliaSmoothOptimizers.github.io/StructuredQR.jl",
    assets = ["assets/style.css"],
  ),
  pages = ["Home" => "index.md", "Reference" => "reference.md"],
)

deploydocs(;
  repo = "github.com/JuliaSmoothOptimizers/StructuredQR.jl",
  push_preview = true,
  devbranch = "main",
)
