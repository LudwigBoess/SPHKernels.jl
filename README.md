| **Documentation**                                                 | **Build Status**                                                                                | **Licence**                                                                                | **Version Status** |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:|:-----------:|
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/SPHKernels.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/SPHKernels.jl/dev) | [![Run CI on master](https://github.com/LudwigBoess/SPHKernels.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml/badge.svg)](https://github.com/LudwigBoess/SPHKernels.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml) [![codecov.io](https://codecov.io/gh/LudwigBoess/SPHKernels.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/SPHKernels.jl?branch=master) | [![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md) | ![TagBot](https://github.com/LudwigBoess/SPHKernels.jl/workflows/TagBot/badge.svg) ![CompatHelper](https://github.com/LudwigBoess/SPHKernels.jl/workflows/CompatHelper/badge.svg) |

# SPHKernels.jl

This package supplies a number of kernels frequently used in Smoothed-Particle Hydrodynamics (SPH), as well as functions to evaluate their values and derivatives in 2D and 3D.

These kernels include the B-splines (`Cubic` and `Quintic`) suggested in [Monaghan & Lattanzio (1985)](https://ui.adsabs.harvard.edu/abs/1985A%26A...149..135M/abstract) and the Wendland functions (`WendlandC2`, `WendlandC4`, `WendlandC6` and `WendlandC8` from [Wendland (2009)](https://www.researchgate.net/publication/220179293_Divergence-Free_Kernel_Methods_for_Approximating_the_Stokes_Problem)) as suggested in [Dehnen & Aly (2012)](https://academic.oup.com/mnras/article/425/2/1068/1187211).

In this implementation we follow the convention of Dehnen&Aly in using the 'compact kernel support' as a means to define the maximum extent of the kernel. They denote this ``H`` in their paper, for convenience (aka for not having to type caps) we use `h` in the code.

## Evaluating Kernels

To evaluate a 3D kernel you need to use the function

```julia
kernel_value(k::SPHKernel, u::Real, h_inv::Real)
```

where [SPHKernel](@ref) is the supertype for an implemented SPH kernel, ``u = \frac{x}{h}`` is the distance to the kernel origin in measures of the compact kernel support and `h_inv` is the inverse of the compact kernel support.

The same goes for a 1D kernel

```julia
kernel_value(k::SPHKernel, u::Real, h_inv::Real)
```


If you want your code to look a little more fancy you can also use the alternative functions [𝒲₁](@ref), where the respective subscript refers to the dimension:

```julia
𝒲( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value(kernel, u, h_inv)
```

As an example:
```julia
using SPHKernels 

# Wendland C6 kernel with double precision in 3D
k     = WendlandC6(Float64, 3)
# distance between the particle and the origin of the kernel
r     = 0.5
h     = 1.0
h_inv = 1.0/h
u     = r * h_inv

# kernel value at position r
val = 𝒲(k, u, h_inv)
```


## Evaluating Derivatives

Similar to [Evaluating Kernels](@ref) you can evluate a kernel derivative with

```julia
kernel_deriv(k::SPHKernel, u::Real, h_inv::Real)
```

or in the fancy way:

```julia
d𝒲(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv(kernel, u, h_inv)
```

## Bias Correction

You can correct for the kernel bias of the Wendland kernels as described in [Dehnen & Aly (2012)](https://academic.oup.com/mnras/article/425/2/1068/1187211), Eq. 18 + 19 with the functions:

```julia
bias_correction(kernel::SPHKernel, density::Real, m::Real, h_inv::Real, n_neighbours::Integer)
```

or again in the fancy way

```julia
δρ(kernel::SPHKernel, density::Real, m::Real, h_inv::Real, n_neighbours::Integer) = bias_correction(kernel, density, m, h_inv, n_neighbours)

```

This will return a new value for the density:

```julia
using SPHKernels
density = 1.0
kernel  = WendlandC6(3)

# correct density
density = bias_correction_3D(kernel, density, 1.0, 0.5, 295)
```