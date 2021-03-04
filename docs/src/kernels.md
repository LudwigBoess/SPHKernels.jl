# Kernels

```@meta
CurrentModule = SPHKernels
DocTestSetup = quote
    using SPHKernels
end
```

This package supplies a number of kernels frequently used in Smoothed-Particle Hydrodynamics (SPH), as well as functions to evaluate their values and derivatives in 2D and 3D.

These kernels include the B-splines ([Cubic](@ref) and [Quintic](@ref)) suggested in [Monaghan & Lattanzio (1985)](https://ui.adsabs.harvard.edu/abs/1985A%26A...149..135M/abstract) and the Wendland functions ([WendlandC2](@ref), [WendlandC4](@ref) and [WendlandC6](@ref)) as suggested in [Dehnen & Aly (2012)](https://academic.oup.com/mnras/article/425/2/1068/1187211).

In this implementation we follow the convention of Dehnen&Aly in using the 'compact kernel support' as a means to define the maximum extent of the kernel. They denote this ``H`` in their paper, for convenience (aka for not having to type caps) we use `h` in the code.

## Evaluating Kernels

To evaluate a 3D kernel you need to use the function

```julia
kernel_value_3D(k::SPHKernel, u::Real, h_inv::Real)
```

where [SPHKernel](@ref) is the supertype for an implemented SPH kernel, ``u = \frac{x}{h}`` is the distance to the kernel origin in measures of the compact kernel support and `h_inv` is the inverse of the compact kernel support.

The same goes for a 1D kernel

```julia
kernel_value_1D(k::SPHKernel, u::Real, h_inv::Real)
```

and a 2D kernel

```julia
kernel_value_1D(k::SPHKernel, u::Real, h_inv::Real)
```

If you want your code to look a little more fancy you can also use the alternative functions [𝒲₁](@ref), where the respective subscript refers to the dimension:

```julia
𝒲₁( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value_1D(kernel, u, h_inv)
𝒲₂( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value_2D(kernel, u, h_inv)
𝒲₃( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value_3D(kernel, u, h_inv)
```

As an example:
```@example
using SPHKernels # hide
k     = WendlandC6()
# distance between the particle and the origin of the kernel
r     = 0.5
h     = 1.0
h_inv = 1.0/h
u     = r * h_inv

# kernel value at position r
val = 𝒲₃(k, u, h_inv)

println("val = $val")
```


## Evaluating Derivatives

Similar to [Evaluating Kernels](@ref) you can evluate a kernel derivative with

```julia
kernel_deriv_3D(k::SPHKernel, u::Real, h_inv::Real)
```

or in the fancy way:

```julia
∇𝒲₁(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv_1D(kernel, u, h_inv)
∇𝒲₂(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv_2D(kernel, u, h_inv)
∇𝒲₃(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv_3D(kernel, u, h_inv)

```

## Bias Correction

You can correct for the kernel bias of the Wendland kernels as described in [Dehnen & Aly (2012)](https://academic.oup.com/mnras/article/425/2/1068/1187211), Eq. 18 + 19 with the functions:

```julia
bias_correction_1D(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)
bias_correction_2D(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)
bias_correction_3D(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)
```

or again in the fancy way

```julia
δρ₁(kernel::SPHKernel, density::Real, m::Real, h_inv::Real) = bias_correction_1D(kernel, density, m, h_inv)
δρ₂(kernel::SPHKernel, density::Real, m::Real, h_inv::Real) = bias_correction_2D(kernel, density, m, h_inv)
δρ₃(kernel::SPHKernel, density::Real, m::Real, h_inv::Real) = bias_correction_3D(kernel, density, m, h_inv)

```

This will return a new value for the density:

```@example
using SPHKernels # hide
density = 1.0
kernel  = WendlandC6()

# correct density
density = bias_correction_3D(kernel, density, 1.0, 0.5)

println("density = $density")
```