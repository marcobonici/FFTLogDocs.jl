using Documenter
using FFTLog
using FFTLogDocs
using Plots
using PlotThemes

ENV["GKSwstype"] = "100"

push!(LOAD_PATH,"../src/")

makedocs(
    modules = [FFTLog, FFTLog],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
    sidebar_sitename=true),
    sitename = "FFTLog.jl",
    authors  = "Marco Bonici",
    pages = [
        "Home" => "index.md",
        "Components" => "components.md"
    ]
)

deploydocs(
    repo = "github.com/marcobonici/FFTLogDocs.jl.git",
    devbranch = "main"
)
