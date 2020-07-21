# Kernels

This package supplies a number of kernels frequently used in Smoothed-Particle Hydrodynamics (SPH), as well as functions to evaluate their values and derivatives in 2D and 3D.

These kernels include the B-splines ([Cubic](@ref) and [Quintic](@ref)) suggested in [Monaghan & Lattanzio (1985)](https://ui.adsabs.harvard.edu/abs/1985A%26A...149..135M/abstract) and the Wendland functions ([WendlandC4](@ref) and [WendlandC6](@ref)) as suggested in [Dehnen & Aly (2012)](https://academic.oup.com/mnras/article/425/2/1068/1187211).

In this implementation we follow the convention of Dehnen&Aly in using the 'compact kernel support' as a means to define the maximum extent of the kernel. They denote this ``H`` in their paper, for convenience (aka for not having to type caps) we use `h` in the code.

## Evaluating Kernels

To evaluate a 2D kernel you need to use the function

```julia
kernel_value_2D(k::SPHKernel, u::Real, h_inv::Real)
```

where [SPHKernel](@ref) is the supertype for an implemented SPH kernel, ``u = \frac{x}{h}`` is the distance to the kernel origin in measures of the compact kernel support and `h_inv` is the inverse of the compact kernel support.

The same goes for a 3D kernel

```julia
kernel_value_3D(k::SPHKernel, u::Real, h_inv::Real)
```

## Evaluating Derivatives

Similar to [Evaluating Kernels](@ref) you can evluate a kernel derivative with

```julia
kernel_deriv_2D(k::SPHKernel, u::Real, h_inv::Real)
```