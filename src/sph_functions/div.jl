"""
    kernel_div( k::AbstractSPHKernel, h_inv::T1, 
                xᵢ::T2, xⱼ::T2, Aⱼ::T2) where {T1,T2}

Compute the kernel divergence `∇⋅𝒲` between particle `i` and neighbour `j` for some SPH quantity `A`.
"""
function kernel_div( k::AbstractSPHKernel, h_inv::T1, 
                     xᵢ::T2, xⱼ::T2, Aⱼ::T2) where {T1,T2}
    Aⱼ ⋅ ∇𝒲(k, h_inv, xᵢ, xⱼ)
end


"""
    ∇dot𝒲( k::AbstractSPHKernel, h_inv::T1, 
            xᵢ::T2, xⱼ::T2, Aⱼ::T2) where {T1,T2}

Compute the kernel divergence `∇⋅𝒲` between particle `i` and neighbour `j` for some SPH quantity `A`.
Compact notation of [`kernel_div`](@ref).
"""
∇dot𝒲( k::AbstractSPHKernel, h_inv::T1, 
     xᵢ::T2, xⱼ::T2, Aⱼ::T2) where {T1,T2} = kernel_div( k, h_inv, xᵢ, xⱼ, Aⱼ)


"""
    quantity_divergence( k::AbstractSPHKernel, h_inv::Real, 
                          xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                          Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
                          mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the divergence of the SPH quantity `A` for particle `i`.

``∇\\cdot\\vec{A}_i(x) ≈ \\sum_j \\frac{m_j}{\\rho_j} \\vec{A}_j \\cdot ∇W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function quantity_divergence( k::AbstractSPHKernel, h_inv::Real, 
                              xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                              Aⱼ::Vector{<:Real},
                              mⱼ::Real,             ρⱼ::Real ) 
     
    mⱼ / ρⱼ * ∇dot𝒲( k, h_inv, xᵢ, xⱼ, Aⱼ)

end

"""
    ∇dotA( k::AbstractSPHKernel, h_inv::Real, 
           xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
           Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
           mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the divergence of the SPH quantity `A` for particle `i`.
Compact notation of [`quantity_divergence`](@ref).

``∇\\cdot\\vec{A}_i(x) ≈ \\sum_j \\frac{m_j}{\\rho_j} \\vec{A}_j \\cdot ∇W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
∇dot𝒜( k::AbstractSPHKernel, h_inv::Real, 
        xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
        Aⱼ::Vector{<:Real},
        mⱼ::Real,             ρⱼ::Real ) = quantity_divergence( k, h_inv, xᵢ, xⱼ, Aⱼ, mⱼ, ρⱼ )

