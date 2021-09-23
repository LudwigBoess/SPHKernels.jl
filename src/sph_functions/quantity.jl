"""
    kernel_value( k::AbstractSPHKernel, h_inv::Real, 
                  xᵢ::Union{Real, Vector{<:Real}}, 
                  xⱼ::Union{Real, Vector{<:Real}} )

Computes the value of the kernel `k` at the position of the neighbour `xⱼ`. 

``W(\vec{x}_i - \vec{x}_j, h_i)``
"""
function kernel_value( k::AbstractSPHKernel, h_inv::Real, 
                       xᵢ::Union{Real, Vector{<:Real}}, 
                       xⱼ::Union{Real, Vector{<:Real}} ) 
    
    r = √( sum( (xᵢ .- xⱼ).^2 ) )

    q  = r*h_inv
    Wq = 𝒲(k, q, h_inv)

    A = map(xᵢ, xⱼ) do i, j
        d = i - j
        Wq*(d/r)
    end
    
    return A
end

"""
    𝒜( k::AbstractSPHKernel, h_inv::Real, 
        xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
        Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
        mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the SPH quantity `A` for particle `i`.

See e.g. Price 2012:
``\vec{A}_i(x) ≈ \sum_j m_j \frac{\vec{A}_j}{\rho_j} W(\vec{x}_i - \vec{x}_j, h_i)``
"""
𝒜( k::AbstractSPHKernel, h_inv::Real, 
    xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
    Aⱼ::Vector{<:Real},
    mⱼ::Real,             ρⱼ::Real ) = mⱼ / ρⱼ * Aⱼ * kernel_value( k, h_inv, xᵢ, xⱼ)