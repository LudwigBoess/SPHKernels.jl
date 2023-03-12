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


    """
        Precompile Functions
    """

    using SnoopPrecompile    # this is a small dependency

    @precompile_setup begin
        # Putting some things in `setup` can reduce the size of the
        # precompile file and potentially make loading faster.
        kernels = [kernel(dt, dim) for kernel ∈ [Cubic, Quintic, WendlandC2, WendlandC4, WendlandC6, WendlandC8, DoubleCosine], dt ∈ [Float32, Float64], dim ∈ [1, 2, 3]]

        @precompile_all_calls begin
            # all calls in this block will be precompiled, regardless of whether
            # they belong to your package or not (on Julia 1.8 and higher)

            for k ∈ kernels
                for u ∈ [ 0.3, 0.5, 0.8, 1.5]
                    kernel_value(k, u, 0.5)
                    𝒲(k, u, 0.5)
                    kernel_deriv(k, u, 0.5)
                    d𝒲(k, u, 0.5)
                    bias_correction(k, 1.0, 1.0, 0.5, 128)
                    δρ(k, 1.0, 1.0, 0.5, 128)
                end
            end

            # test quantities setup
            x_i = [0.0, 0.0, 0.0]
            x_j = [0.5, 0.5, 0.5]
            Δx = x_i - x_j
            A_i = [1.0, 1.0, 1.0]
            A_j = [1.5, 1.5, 1.5]
            m_j = 1.5
            ρ_j = 1.5
            k = WendlandC6(3)
            r = SPHKernels.get_r(x_i, x_j)
            h_inv = 1.0
            u = r * h_inv

            for k ∈ kernels
                # quantity
                𝒲(k, h_inv, x_i, x_j)
                𝒲(k, u, h_inv)
                𝒜(k, h_inv, x_i, x_j, A_j[1], m_j, ρ_j)
                𝒜(k, r, h_inv, A_j[1], m_j, ρ_j)

                # gradient
                ∇𝒲(k, h_inv, x_i, x_j)
                ∇𝒲(k, h_inv, x_i[1], x_j[1])
                ∇𝒲(k, r, h_inv, Δx)
                ∇𝒜(k, h_inv, x_i, x_j, A_j, m_j, ρ_j)
                ∇𝒜(k, h_inv, x_i[1], x_j[1], A_j[1], m_j, ρ_j)
                ∇𝒜(k, r, h_inv, Δx, A_j, m_j, ρ_j)

                # divergence
                ∇̇dot𝒲(k, h_inv, x_i, x_j, A_j)
                ∇dot𝒜(k, h_inv, x_i, x_j, A_j, m_j, ρ_j)

                # curl 
                ∇x𝒲(k, h_inv, x_i, x_j, A_j)
                ∇x𝒜(k, h_inv, x_i, x_j, A_j, m_j, ρ_j)
            end


        end
    end
    
end # module
