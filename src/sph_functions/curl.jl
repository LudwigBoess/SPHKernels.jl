"""
    kernel_curl(k::AbstractSPHKernel, h_inv::T1, xᵢ::T2, xⱼ::T2, Aⱼ::T2) where {T1,T2}

Compute the kernel curl `∇x𝒲` between particle `i` and neighbour `j` for some SPH quantity `A`.
"""
function kernel_curl( k::AbstractSPHKernel, h_inv::T1, 
                      xᵢ::T2, xⱼ::T2, Aⱼ::T2) where {T1,T2}
    Aⱼ × ∇𝒲(k, h_inv, xᵢ, xⱼ)
end


"""
    ∇x𝒲(k::AbstractSPHKernel, h_inv::T1, xᵢ::T2, xⱼ::T2, Aⱼ::T2) where {T1,T2}

Compute the kernel curl `∇x𝒲` between particle `i` and neighbour `j` for some SPH quantity `A`.
"""
∇x𝒲(k::AbstractSPHKernel, h_inv::T1, xᵢ::T2, xⱼ::T2, Aⱼ::T2) where {T1,T2} = kernel_curl( k, h_inv, xᵢ, xⱼ, Aⱼ)

"""
    quantity_curl(k::AbstractSPHKernel, h_inv::T1, xᵢ::T2, xⱼ::T2, Aⱼ::T2, mⱼ::T1, ρⱼ::T1 ) where {T1,T2}

Compute the contribution of particle `j` to the curl of the SPH quantity `A` for particle `i`.

``∇×\\vec{A}_i(x) ≈ - \\sum_j \\frac{m_j}{\\rho_j} \\vec{A}_j \\times ∇W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function quantity_curl( k::AbstractSPHKernel, h_inv::T1, 
                        xᵢ::T2, xⱼ::T2, Aⱼ::T2,
                        mⱼ::T1, ρⱼ::T1 ) where {T1,T2}
    (-mⱼ / ρⱼ) .* ∇x𝒲( k, h_inv, xᵢ, xⱼ, Aⱼ)
end


"""
    ∇x𝒜(k::AbstractSPHKernel, h_inv::T1, xᵢ::T2, xⱼ::T2, Aⱼ::T2, mⱼ::T1, ρⱼ::T1 ) where {T1,T2}

Compute the contribution of particle `j` to the curl of the SPH quantity `A` for particle `i`.

``∇×\\vec{A}_i(x) ≈ - \\sum_j \\frac{m_j}{\\rho_j} \\vec{A}_j \\times ∇W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
∇x𝒜(k::AbstractSPHKernel, h_inv::T1, xᵢ::T2, xⱼ::T2, Aⱼ::T2, mⱼ::T1, ρⱼ::T1 ) where {T1,T2} = 
    quantity_curl( k, h_inv, xᵢ, xⱼ, Aⱼ, mⱼ, ρⱼ)