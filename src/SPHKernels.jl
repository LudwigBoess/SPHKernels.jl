"""
    This package contains a number of kernels frequently used in SPH, as well as their functions and derivatives an 2D and 3D.
"""

module SPHKernels

    export  kernel_value,     𝒲,
            kernel_deriv,    d𝒲,
            kernel_gradient, ∇𝒲,
            bias_correction, δρ,
            AbstractSPHKernel,
            Cubic, 
            Quintic,
            WendlandC2,
            WendlandC4,
            WendlandC6,
            WendlandC8
        

    """
        AbstractSPHKernel

    Supertype for all SPH kernels.
    """
    abstract type AbstractSPHKernel end

    include("bsplines/Cubic.jl")
    include("bsplines/Quintic.jl")
    include("wendland/C2.jl")
    include("wendland/C4.jl")
    include("wendland/C6.jl")
    include("wendland/C8.jl")
    include("sph_functions/gradients.jl")

    # multiple dispatch for nicer look

    """
        𝒲₁( kernel::AbstractSPHKernel, u::Real, h_inv::Real)

    Evaluate 1D spline at position ``u = \\frac{x}{h}``.
    """
    𝒲( kernel::AbstractSPHKernel, u::Real, h_inv::Real) = kernel_value(kernel, u, h_inv)


    """
        d𝒲₂( kernel::AbstractSPHKernel, u::Real, h_inv::Real)

    Evaluate 1D derivative at position ``u = \\frac{x}{h}``.
    """
    d𝒲(kernel::AbstractSPHKernel, u::Real, h_inv::Real) = kernel_deriv(kernel, u, h_inv)

    """ 
        δρ₁(kernel::AbstractSPHKernel, density::Real, m::Real, h_inv::Real)

    Corrects the 1D density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
    """
    δρ(kernel::AbstractSPHKernel, density::Real, m::Real, h_inv::Real, n_neighbours::Integer) = bias_correction(kernel, density, m, h_inv, n_neighbours)


    
end # module
