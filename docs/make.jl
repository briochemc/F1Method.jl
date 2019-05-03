using Documenter, F1Method

makedocs(
    sitename="F1Method Documentation",
    # options
    modules = [F1Method]
)

deploydocs(
    repo = "github.com/briochemc/F1Method.jl.git",
)