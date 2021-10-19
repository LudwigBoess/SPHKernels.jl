"""
    kernel_value( k::AbstractSPHKernel, h_inv::Real, 
                  xᵢ::Real, xⱼ::Real )

Computes the value of the kernel `k` at the position of the neighbour `xⱼ`. 

``W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function kernel_value( k::AbstractSPHKernel, h_inv::Real, 
                       xᵢ::Real, xⱼ::Real ) 
    
    u  = abs(xᵢ - xⱼ)*h_inv

    𝒲(k, u, h_inv)
end

"""
    kernel_value( k::AbstractSPHKernel, h_inv::Real, 
                  xᵢ::Union{Real, Vector{<:Real}}, 
                  xⱼ::Union{Real, Vector{<:Real}} )

Computes the value of the kernel `k` at the position of the neighbour `xⱼ`. 

``W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function kernel_value( k::AbstractSPHKernel, h_inv::Real, 
                       xᵢ::Vector{<:Real}, xⱼ::Vector{<:Real} )
    
    u  = get_r(xᵢ, xⱼ) * h_inv

    𝒲(k, u, h_inv)
end


"""
    𝒜( k::AbstractSPHKernel, h_inv::Real, 
        xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
        Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
        mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the SPH quantity `A` for particle `i`.

See e.g. Price 2012:
``\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
𝒲( k::AbstractSPHKernel, h_inv::Real, 
    xᵢ::Union{Real, Vector{<:Real}}, 
    xⱼ::Union{Real, Vector{<:Real}}) = kernel_value( k, h_inv, xᵢ, xⱼ)

"""
    kernel_quantity( k::AbstractSPHKernel, h_inv::Real, 
                     xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                     Aⱼ::Vector{<:Real},
                     mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the SPH quantity `A` for particle `i`.
Based on Euclidean distance `r` between the particles. 
Useful if many quantities need to be computed for the same particle pair.

``\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function kernel_quantity( k::AbstractSPHKernel, 
                          r::Real, h_inv::Real, 
                          Aⱼ::Union{Real, Vector{<:Real}},
                          mⱼ::Real,             ρⱼ::Real ) 
                 
    mj_wk = mⱼ / ρⱼ * 𝒲(k, r*h_inv, h_inv)

    map(Aⱼ) do Aj
        mj_wk * Aj
    end
end


"""
    kernel_quantity( k::AbstractSPHKernel, h_inv::Real, 
                     xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                     Aⱼ::Vector{<:Real},
                     mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the SPH quantity `A` for particle `i`.
Based on positions `xᵢ` and `xⱼ`.

``\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function kernel_quantity( k::AbstractSPHKernel, h_inv::Real, 
                 xᵢ::Union{Real, Vector{<:Real}}, 
                 xⱼ::Union{Real, Vector{<:Real}},
                 Aⱼ::Union{Real, Vector{<:Real}},
                 mⱼ::Real,             ρⱼ::Real )
                 
    mj_wk = mⱼ / ρⱼ * kernel_value( k, h_inv, xᵢ, xⱼ)

    map(Aⱼ) do Aj
        mj_wk * Aj
    end
end


"""
    𝒜( k::AbstractSPHKernel, h_inv::Real, 
        xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
        Aⱼ::Vector{<:Real},
        mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the SPH quantity `A` for particle `i`.
Based on positions `xᵢ` and `xⱼ`.

``\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
𝒜( k::AbstractSPHKernel, h_inv::Real, 
    xᵢ::Union{Real, Vector{<:Real}}, 
    xⱼ::Union{Real, Vector{<:Real}},
    Aⱼ::Union{Real, Vector{<:Real}},
    mⱼ::Real,             ρⱼ::Real ) = kernel_quantity( k, h_inv, xᵢ, xⱼ, Aⱼ, mⱼ, ρⱼ)


"""
    𝒜( k::AbstractSPHKernel, 
        r::Real,  h_inv::Real, 
        Aⱼ::Union{Real, Vector{<:Real}},
        mⱼ::Real, ρⱼ::Real )

Compute the contribution of particle `j` to the SPH quantity `A` for particle `i`.
Based on Euclidean distance `r` between the particles. 
Useful if many quantities need to be computed for the same particle pair.

``\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} W(r, h_i)``
"""
𝒜( k::AbstractSPHKernel, 
    r::Real,  h_inv::Real, 
    Aⱼ::Union{Real, Vector{<:Real}},
    mⱼ::Real, ρⱼ::Real ) = kernel_quantity( k, r, h_inv, Aⱼ, mⱼ, ρⱼ)