"""
    kernel_div( k::SPHKernel,       h_inv::Real, 
                xᵢ::Vector{<:Real}, xⱼ::Vector{<:Real},
                Aᵢ::Vector{<:Real}, Aⱼ::Vector{<:Real} )

Compute the kernel divergence `∇⋅𝒲` between particle `i` and neighbour `j` for some SPH quantity `A`.
"""
function kernel_div( k::SPHKernel,       h_inv::Real, 
                     xᵢ::Vector{<:Real}, xⱼ::Vector{<:Real},
                     Aᵢ::Vector{<:Real}, Aⱼ::Vector{<:Real})
    
    (Aᵢ - Aⱼ) ⋅ ∇𝒲(k, h_inv, xᵢ, xⱼ)
end

"""
    ∇̇dot𝒲( k::AbstractSPHKernel, h_inv::Real, 
            xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
            Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real})

Compute the kernel divergence `∇⋅𝒲` between particle `i` and neighbour `j` for some SPH quantity `A`.
Compact notation of [`kernel_div`](@ref).
"""
∇̇dot𝒲( k::AbstractSPHKernel, h_inv::Real, 
     xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
     Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real}) = kernel_div( k, h_inv, 
                                                              xᵢ, xⱼ, 
                                                              Aᵢ, Aⱼ)

"""
    quantity_divergence( k::AbstractSPHKernel, h_inv::Real, 
                          xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                          Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
                          mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the divergence of the SPH quantity `A` for particle `i`.

``∇\cdot\vec{A}_i(x) ≈ \sum_j m_j \frac{\vec{A}_j}{\rho_j} \cdot ∇W(\vec{x}_i - \vec{x}_j, h_i)``
"""
function quantity_divergence( k::AbstractSPHKernel, h_inv::Real, 
                               xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                               Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
                               mⱼ::Real,             ρⱼ::Real ) 
     
    mⱼ / ρⱼ * ∇̇dot𝒲( k, h_inv, xᵢ, xⱼ, Aᵢ, Aⱼ)

end

"""
    ∇dotA( k::AbstractSPHKernel, h_inv::Real, 
           xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
           Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
           mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the divergence of the SPH quantity `A` for particle `i`.
Compact notation of [`quantity_divergence`](@ref).

``∇\cdot\vec{A}_i(x) ≈ \sum_j m_j \frac{\vec{A}_j}{\rho_j} \cdot ∇W(\vec{x}_i - \vec{x}_j, h_i)``
"""
∇dot𝒜( k::AbstractSPHKernel, h_inv::Real, 
        xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
        Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
        mⱼ::Real,             ρⱼ::Real ) = quantity_divergence( k, h_inv, xᵢ, xⱼ, Aᵢ, Aⱼ, mⱼ, ρⱼ )