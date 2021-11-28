using Documenter
using SPHKernels

makedocs(
    sitename = "SPHKernels",
    format = Documenter.HTML(),
    modules = [SPHKernels],
    pages = [
            "Table of Contents" => "index.md",
            "Install" => "install.md",
            "Kernels" => "kernels.md",
            "Extending" => "extending.md",
            "SPH Functions" => "sph.md",
            "API reference" => "api.md"
            ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/LudwigBoess/SPHKernels.jl.git"
)
