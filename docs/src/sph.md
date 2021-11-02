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

## Quantities

You can compute the contribution of particle `j` to the SPH quantity `A` for particle `i` via (see e.g. Price 2012):

``\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} W(\\vec{x}_i - \\vec{x}_j, h_i)``

To do this you can loop over the function

```@docs
kernel_quantity
```

or its more compact form

```@docs
𝒜
```

## Gradient

The gradient of a quantity in SPH can be reduced to the gradient of the kernel (see e.g. Price 2012):

``∇\vec{A}_i(x) ≈ \sum_j m_j \frac{\vec{A}_j}{\rho_j} ∇W(||\vec{x}_i - \vec{x}_j||, h_i)``

We provide two functionalities 

### Kernel

You can compute the gradient of the kernel at position `x_j` 
``∇W(x_{ij}, h_i) = \frac{dW}{dx}\vert_{x_j} \frac{Δx_{ij}}{||x_{ij}||} \frac{1}{h_i}``

by using

```@docs
kernel_gradient
```

or its more compact form

```@docs
∇𝒲
```

where `xᵢ` and `xⱼ` are the positions of particles `i` and `j` in 1D-3D space.

### Quantity

To compute the gradient of a SPH quantity for particle `i` 
``∇\vec{A}_i(x) ≈ - \sum_j m_j \frac{\vec{A}_j}{\rho_j} ∇W(\vec{x}_i - \vec{x}_j, h_i)``

you can loop over

```@docs
quantity_gradient
```

or its compact form

```@docs
∇𝒜
```

## Divergence

The divergence of a quantity in SPH can be calculated as (see e.g. Price 2012):

``∇\\cdot\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} \\cdot ∇W(\\vec{x}_i - \\vec{x}_j, h_i)``

We provide two functionalities 

### Kernel

You can compute the divergence of the kernel at position `x_j` by using

```@docs
kernel_div
```

or its more compact form

```@docs
∇̇dot𝒲
```

where `xᵢ` and `xⱼ` are the positions of particles `i` and `j` in 1D-3D space.

### Quantity

To compute the gradient of a SPH quantity for particle `i` you can loop over

```@docs
quantity_divergence
```

or its compact form

```@docs
∇dot𝒜
```

## Curl

The curl of a quantity in SPH can be calculated as (see e.g. Price 2012):

``∇×\\vec{A}_i(x) ≈ - \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} \\times ∇W(\\vec{x}_i - \\vec{x}_j, h_i)``

We provide two functionalities 

### Kernel

You can compute the divergence of the kernel at position `x_j` by using

```@docs
kernel_curl
```

or its more compact form

```@docs
∇x𝒲
```

where `xᵢ` and `xⱼ` are the positions of particles `i` and `j` in 1D-3D space.

### Quantity

To compute the gradient of a SPH quantity for particle `i` you can loop over

```@docs
quantity_curl
```

or its compact form

```@docs
∇x𝒜
```