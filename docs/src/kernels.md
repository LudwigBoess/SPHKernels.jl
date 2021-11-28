# Kernels

```@meta
CurrentModule = SPHKernels
DocTestSetup = quote
    using SPHKernels
end
```

This package supplies a number of kernels frequently used in Smoothed-Particle Hydrodynamics (SPH), as well as functions to evaluate their values and derivatives in 2D and 3D.

These kernels include the B-splines ([Cubic](@ref) and [Quintic](@ref)) suggested in [Monaghan & Lattanzio (1985)](https://ui.adsabs.harvard.edu/abs/1985A%26A...149..135M/abstract) and the Wendland functions ([WendlandC2](@ref), [WendlandC4](@ref), [WendlandC6](@ref)) and [WendlandC8](@ref) ([Wendland 2009](https://www.researchgate.net/publication/220179293_Divergence-Free_Kernel_Methods_for_Approximating_the_Stokes_Problem)) as suggested in [Dehnen & Aly (2012)](https://academic.oup.com/mnras/article/425/2/1068/1187211).

In this implementation we follow the convention of Dehnen&Aly in using the 'compact kernel support' as a means to define the maximum extent of the kernel. They denote this ``H`` in their paper, for convenience (aka for not having to type caps) we use `h` in the code.

```@eval 
using CairoMakie
import ColorSchemes

function get_kernel_values(x, k)            

    W = Vector{Float64}(undef, length(x))
    for i = 1:length(x)
        W[i] = 𝒲(k, x[i], 1.0)
    end

    W
end

function get_kernel_deriv(x, k)            

    dW = Vector{Float64}(undef, length(x))
    for i = 1:length(x)
        dW[i] = d𝒲(k, x[i], 1.0)
    end

    dW
end


function get_kernels(dim)

    [ Cubic(dim), Quintic(dim), 
      WendlandC2(dim), WendlandC4(dim), WendlandC6(dim), WendlandC8(dim)
    ]

end

N_samples = 1_000
x = LinRange(0, 1, N_samples)

dims = [1, 2, 3]

labels = ["Cubic", "Quintic", 
         "Wendland C2", "Wendland C4", "Wendland C6", "Wendland C8"]

colors = [ColorSchemes.BuPu_7[end], ColorSchemes.BuPu_7[end-1], 
          ColorSchemes.PuBuGn_9[end], ColorSchemes.PuBuGn_9[end-1], ColorSchemes.PuBuGn_9[end-2], ColorSchemes.PuBuGn_9[end-3]]

fs    = 25
scale = 750
fig   = Figure(resolution = (1.5*scale, 2.1*scale), fontsize=fs)

for dim = 1:3

    ax_l = Axis(fig[dim+1, 1], xlabel="x = r/h", ylabel="W(x)")
    ax_r = Axis(fig[dim+1, 2], xlabel="x = r/h", ylabel="W'(x)")

    kernels = get_kernels(dim)

    max_W = 0.0

    for i ∈ 1:length(kernels)

        W = get_kernel_values(x, kernels[i]) 
        lines!(ax_l, x, W, label=labels[i], color=colors[i])

        if maximum(W) > max_W 
            max_W = maximum(W)
        end

        dW = get_kernel_deriv(x, kernels[i]) 
        lines!(ax_r, x, dW, label=labels[i], color=colors[i])

    end

    text!(ax_l, "$(dim)D", position = (0.9, 0.95max_W), 
            align = (:center, :center), textsize =2fs)

end

line_elem = [LineElement(color = colors[i], linestyle = nothing) for i = 1:length(colors)]
Legend(fig[1, 1:2], line_elem, labels,
        framevisible = false, 
        orientation = :horizontal,
        nbanks = 2)

save("kernels.png", fig); nothing # hide
```

![kernels](kernels.png)

## Defining Kernels

We use multiple dispatch to make to conform to Julia coding standards and make the code more readable.

To use e.g. a 3D [WendlandC6](@ref) kernel use

```julia
k = WendlandC6(3)
```

This will default to a 3D kernel with a precision defined by the system OS (usually `Float64`).
If you want to use a `Float32` kernel you can define the precision as the (optional) first argument

```julia
k = WendlandC6(Float32, 3)
```

## Evaluating Kernels

To evaluate a kernel you need to use the function

```julia
kernel_value(k::AbstractSPHKernel, u::Real, h_inv::Real)
```

where [AbstractSPHKernel](@ref) is the supertype for an implemented SPH kernel, ``u = \frac{x}{h}`` is the distance to the kernel origin in measures of the compact kernel support and `h_inv` is the inverse of the compact kernel support.

You need to define the dimension of the kernel in the `kernel <: AbstractSPHKernel`, as explained before.

If you want your code to look a little more fancy you can also use the alternative functions [𝒲](@ref).:

```julia
𝒲( kernel::AbstractSPHKernel, u::Real, h_inv::Real) = kernel_value(kernel, u, h_inv)
```

As an example:
```@example
using SPHKernels # hide

# Wendland C6 kernel with double precision in 3D
k     = WendlandC6(Float64, 3)
# distance between the particle and the origin of the kernel
r     = 0.5
h     = 1.0
h_inv = 1.0/h
u     = r * h_inv

# kernel value at position r
val = 𝒲(k, u, h_inv)

println("val = $val")
```


## Evaluating Derivatives

Similar to [Evaluating Kernels](@ref) you can evluate a kernel derivative with

```julia
kernel_deriv(k::AbstractSPHKernel, u::Real, h_inv::Real)
```

or in the fancy way:

```julia
d𝒲(kernel::AbstractSPHKernel, u::Real, h_inv::Real) = kernel_deriv(kernel, u, h_inv)
```

## Bias Correction

You can correct for the kernel bias of the Wendland kernels as described in [Dehnen & Aly (2012)](https://academic.oup.com/mnras/article/425/2/1068/1187211), Eq. 18 + 19 with the functions:

```julia
bias_correction(kernel::AbstractSPHKernel, density::Real, m::Real, h_inv::Real, n_neighbours::Integer)
```

or again in the fancy way

```julia
δρ(kernel::AbstractSPHKernel, density::Real, m::Real, h_inv::Real, n_neighbours::Integer) = bias_correction(kernel, density, m, h_inv, n_neighbours)

```

This will return a new value for the density:

```@example
using SPHKernels # hide
density = 1.0
kernel  = WendlandC6(3)

# correct density
density = bias_correction_3D(kernel, density, 1.0, 0.5, 295)

println("density = $density")
```