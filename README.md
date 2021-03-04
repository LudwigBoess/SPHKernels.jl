| **Documentation**                                                 | **Build Status**                                                                                | **Licence**                                                                                | **Version Status** |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:|:-----------:|
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/SPHKernels.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/SPHKernels.jl/dev) | [![Run CI on master](https://github.com/LudwigBoess/SPHKernels.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml/badge.svg)](https://github.com/LudwigBoess/SPHKernels.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml) [![codecov.io](https://codecov.io/gh/LudwigBoess/SPHKernels.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/SPHKernels.jl?branch=master) | [![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md) | ![TagBot](https://github.com/LudwigBoess/SPHKernels.jl/workflows/TagBot/badge.svg) ![CompatHelper](https://github.com/LudwigBoess/SPHKernels.jl/workflows/CompatHelper/badge.svg) |

# SPHKernels.jl

This package supplies a number of kernels frequently used in Smoothed-Particle Hydrodynamics (SPH), as well as functions to evaluate their values and derivatives in 2D and 3D.

These kernels include the B-splines (`Cubic` and `Quintic`) suggested in [Monaghan & Lattanzio (1985)](https://ui.adsabs.harvard.edu/abs/1985A%26A...149..135M/abstract) and the Wendland functions (`WendlandC2`, `WendlandC4` and `WendlandC6`) as suggested in [Dehnen & Aly (2012)](https://academic.oup.com/mnras/article/425/2/1068/1187211).

In this implementation we follow the convention of Dehnen&Aly in using the 'compact kernel support' as a means to define the maximum extent of the kernel. They denote this ``H`` in their paper, for convenience (aka for not having to type caps) we use `h` in the code.

## Evaluating Kernels

To evaluate a 3D kernel you need to use the function

```julia
kernel_value_3D(k::SPHKernel, u::Real, h_inv::Real)
```

where `SPHKernel` is the supertype for an implemented SPH kernel, ``u = \frac{x}{h}`` is the distance to the kernel origin in measures of the compact kernel support and `h_inv` is the inverse of the compact kernel support.

The same goes for a 1D kernel

```julia
kernel_value_1D(k::SPHKernel, u::Real, h_inv::Real)
```

and a 2D kernel

```julia
kernel_value_2D(k::SPHKernel, u::Real, h_inv::Real)
```

If you want your code to look a little more fancy you can also use the alternative functions `𝒲₁`, where the respective subscript refers to the dimension:

```julia
𝒲₁( kernel::SPHKernel, u::Real, h_inv::Real ) = kernel_value_1D(kernel, u, h_inv)
𝒲₂( kernel::SPHKernel, u::Real, h_inv::Real ) = kernel_value_2D(kernel, u, h_inv)
𝒲₃( kernel::SPHKernel, u::Real, h_inv::Real ) = kernel_value_3D(kernel, u, h_inv)
```

Please see the docs for more details!