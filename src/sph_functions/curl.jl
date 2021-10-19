"""
    kernel_curl( k::AbstractSPHKernel,       h_inv::Real, 
                 xᵢ::Vector{<:Real}, xⱼ::Vector{<:Real},
                 Aᵢ::Vector{<:Real}, Aⱼ::Vector{<:Real} )

Compute the kernel curl `∇x𝒲` between particle `i` and neighbour `j` for some SPH quantity `A`.
"""
function kernel_curl( k::AbstractSPHKernel,       h_inv::Real, 
                      xᵢ::Vector{<:Real}, xⱼ::Vector{<:Real},
                      Aⱼ::Vector{<:Real})
    
    Aⱼ × ∇𝒲(k, h_inv, xᵢ, xⱼ)
end


"""
    ∇x𝒲( k::AbstractSPHKernel, h_inv::Real, 
          xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
          Aⱼ::Vector{<:Real})

Compute the kernel curl `∇x𝒲` between particle `i` and neighbour `j` for some SPH quantity `A`.
"""
∇x𝒲( k::AbstractSPHKernel, h_inv::Real, 
     xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
     Aⱼ::Vector{<:Real}) = kernel_curl( k, h_inv, xᵢ, xⱼ, Aⱼ)

"""
    quantity_curl( k::AbstractSPHKernel, h_inv::Real, 
          xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
          Aⱼ::Vector{<:Real},
          mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the curl of the SPH quantity `A` for particle `i`.

``∇×\\vec{A}_i(x) ≈ - \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} \\times ∇W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function quantity_curl( k::AbstractSPHKernel, h_inv::Real, 
                        xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                        Aⱼ::Vector{<:Real},
                        mⱼ::Real,             ρⱼ::Real ) 
     
    (-mⱼ / ρⱼ) .* ∇x𝒲( k, h_inv, xᵢ, xⱼ, Aⱼ)

end


"""
    ∇x𝒜( k::AbstractSPHKernel, h_inv::Real, 
          xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
          Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
          mⱼ::Real,             ρⱼ::Real )

Compute the contribution of particle `j` to the curl of the SPH quantity `A` for particle `i`.

``∇×\\vec{A}_i(x) ≈ - \\sum_j m_j \\frac{}{\\rho_j} \\times ∇W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
∇x𝒜( k::AbstractSPHKernel, h_inv::Real, 
      xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
      Aⱼ::Vector{<:Real},
      mⱼ::Real,             ρⱼ::Real ) = quantity_curl( k, h_inv, xᵢ, xⱼ, Aⱼ, mⱼ, ρⱼ)