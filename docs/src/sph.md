```@meta
CurrentModule = SPHKernels
DocTestSetup = quote
    using SPHKernels
end
```

# SPH Functions

We provide some basic functions to compute the gradient, divergence and curl of the kernel with a particle quantity if required.

Please note that these functions are not terribly flexible and/or performance optimized, as we believe that this should be done in a seperate `SPHLoop.jl` package. These functions are more meant as a starting point for a more efficient implementation or code that is not used in performance-relevant applications.

## Notation

The derivative of a quantity in gradient, divergence and curl can be reduced to the derivative of the kernel combined with the remaining mathematical operation (see e.g. [Price (2012)](https://ui.adsabs.harvard.edu/abs/2012JCoPh.231..759P/abstract)).

To provide some helper functions we split this into two functions: 

One to calculate only the kernel derivative and use as a gradient for divergence and curl with the particle quantity.

The other to directly calculate the contribution of particle `j` to the gradient, divergence and curl of the SPH quantity for particle `i`.

## Gradient

The gradient of a quantity in SPH can be reduced to the gradient of the kernel (see e.g. Price 2012):

``∇\vec{A}_i(x) ≈ \sum_j m_j \frac{\vec{A}_j}{\rho_j} ∇W(||\vec{x}_i - \vec{x}_j||, h_i)``

We provide two functionalities 

### Kernel

You can compute the gradient of the kernel at position `x_j` 
``∇W(x_{ij}, h_i) = \frac{dW}{dx}\vert_{x_j} \frac{Δx_{ij}}{||x_{ij}||} \frac{1}{h_i}`` 
by using

```julia
kernel_gradient( k::SPHKernel, h_inv::Real, 
                 xᵢ::Union{Real, Vector{<:Real}}, 
                 xⱼ::Union{Real, Vector{<:Real}} )
```

or its more compact form

```julia
∇𝒲( k::SPHKernel, h_inv::Real, 
     xᵢ::Union{Real, Vector{<:Real}}, 
     xⱼ::Union{Real, Vector{<:Real}} )
```
where `xᵢ` and `xⱼ` are the positions of particles `i` and `j` in 1D-3D space.

### Quantity

To compute the gradient of a SPH quantity for particle `i` 
``∇\vec{A}_i(x) ≈ - \sum_j m_j \frac{\vec{A}_j}{\rho_j} ∇W(\vec{x}_i - \vec{x}_j, h_i)``
you can loop over 

```julia
quantity_gradient( k::AbstractSPHKernel, h_inv::Real, 
                   xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                   Aⱼ::Vector{<:Real},
                   mⱼ::Real,             ρⱼ::Real )
```

or its compact form

```julia
∇𝒜( k::AbstractSPHKernel, h_inv::Real, 
     xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
     Aⱼ::Vector{<:Real},
     mⱼ::Real,             ρⱼ::Real )
```

## Divergence