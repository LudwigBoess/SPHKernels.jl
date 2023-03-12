"""
    This package contains a number of kernels frequently used in SPH, as well as their functions and derivatives an 2D and 3D.
"""

module SPHKernels

    export  kernel_value,          𝒲,  𝒜,
            kernel_deriv,         d𝒲,
            bias_correction,       δρ,
            kernel_gradient,      ∇𝒲, 
            quantity_gradient,    ∇𝒜,
            kernel_div,           ∇̇dot𝒲,
            quantity_divergence,  ∇dot𝒜,
            kernel_curl,          ∇x𝒲,
            quantity_curl,        ∇x𝒜,
            AbstractSPHKernel, 
            Cubic, 
            Quintic,
            WendlandC2,
            WendlandC4,
            WendlandC6,
            WendlandC8,
            Tophat,
            DoubleCosine
            
    using LinearAlgebra

    """
        AbstractSPHKernel

    Supertype for all SPH kernels.
    """
    abstract type AbstractSPHKernel end

    """
        get_r(xᵢ::Vector{<:Real}, xⱼ::Vector{<:Real})

    Eukledian distance between `xᵢ` and `xⱼ`.
    """
    function get_r(xᵢ::Vector{<:Real}, xⱼ::Vector{<:Real})
        # eukledian distance
        r2 = 0
        @inbounds for dim ∈ eachindex(xᵢ)
            r2 += (xᵢ[dim] - xⱼ[dim])^2
        end
        √(r2)
    end

    """
        get_r(xᵢ::Real, xⱼ::Real)

    Eukledian distance between `xᵢ` and `xⱼ`.
    """
    get_r(xᵢ::Real, xⱼ::Real) = abs(xᵢ - xⱼ)

    include("bsplines/Cubic.jl")
    include("bsplines/Quintic.jl")
    include("wendland/C2.jl")
    include("wendland/C4.jl")
    include("wendland/C6.jl")
    include("wendland/C8.jl")
    include("tophat/tophat.jl")
    include("trigonometric/double_cosine.jl")
    include("sph_functions/gradient.jl")
    include("sph_functions/div.jl")
    include("sph_functions/curl.jl")
    include("sph_functions/quantity.jl")

    # multiple dispatch for nicer look
    """
        𝒲( kernel::AbstractSPHKernel, u::Real, h_inv::Real)

    Evaluate kernel at position ``u = \\frac{x}{h}``.
    """
    𝒲( kernel::AbstractSPHKernel, u::Real, h_inv::Real) = kernel_value(kernel, u, h_inv)


    """
        d𝒲( kernel::AbstractSPHKernel, u::Real, h_inv::Real)

    Evaluate derivative at position ``u = \\frac{x}{h}``.
    """
    d𝒲(kernel::AbstractSPHKernel, u::Real, h_inv::Real) = kernel_deriv(kernel, u, h_inv)

    """ 
        δρ₁(kernel::AbstractSPHKernel, density::Real, m::Real, h_inv::Real)

    Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
    """
    δρ(kernel::AbstractSPHKernel, density::Real, m::Real, h_inv::Real, n_neighbours::Integer) = bias_correction(kernel, density, m, h_inv, n_neighbours)

    
end # module
