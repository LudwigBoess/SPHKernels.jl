"""
    kernel_value( k::AbstractSPHKernel, h_inv::Real, 
                  xᵢ::Real, xⱼ::Real )

Computes the value of the kernel `k` at the position of the neighbour `xⱼ`. 

``W(x_i - x_j, h_i)``
"""
function kernel_value( k::AbstractSPHKernel, h_inv::Real, 
                       xᵢ::Real, xⱼ::Real ) 
    
    u  = abs(xᵢ - xⱼ)*h_inv

    𝒲(k, u, h_inv)
end

"""
    kernel_value( k::AbstractSPHKernel, h_inv::T1, 
                       xᵢ::T2, xⱼ::T2 ) where {T1,T2}

Computes the value of the kernel `k` at the position of the neighbour `xⱼ`. 

``W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function kernel_value( k::AbstractSPHKernel, h_inv::T1, 
                       xᵢ::T2, xⱼ::T2 ) where {T1,T2}
    
    u  = get_r(xᵢ, xⱼ) * h_inv

    𝒲(k, u, h_inv)
end


"""
    𝒲( k::AbstractSPHKernel, h_inv, xᵢ, xⱼ)

Computes the value of the kernel `k` at the position of the neighbour `xⱼ`. 

``W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
𝒲( k::AbstractSPHKernel, h_inv, xᵢ, xⱼ) = kernel_value( k, h_inv, xᵢ, xⱼ)

"""
    kernel_quantity(k::AbstractSPHKernel, r::T1, h_inv::T1, 
                    Aⱼ::T2, mⱼ::T1, ρⱼ::T1 ) where {T1,T2}

Compute the contribution of particle `j` to the SPH quantity `A` for particle `i`.
Based on Euclidean distance `r` between the particles. 
Useful if many quantities need to be computed for the same particle pair.

``\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function kernel_quantity( k::AbstractSPHKernel, r::T1, h_inv::T1, 
                          Aⱼ::T2, mⱼ::T1, ρⱼ::T1 ) where {T1,T2}
                 
    mj_wk = mⱼ / ρⱼ * 𝒲(k, r*h_inv, h_inv)

    map(Aⱼ) do Aj
        mj_wk * Aj
    end
end


"""
    kernel_quantity(k::AbstractSPHKernel, h_inv::T1, 
                    xᵢ::T2, xⱼ::T2, Aⱼ::T2, mⱼ::T1, ρⱼ::T1 ) where {T1,T2}

Compute the contribution of particle `j` to the SPH quantity `A` for particle `i`.
Based on positions `xᵢ` and `xⱼ`.

``\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
function kernel_quantity(k::AbstractSPHKernel, h_inv::T1, 
                         xᵢ::T2, xⱼ::T2, Aⱼ::Union{T1,T2}, mⱼ::T1, ρⱼ::T1 ) where {T1,T2}
                 
    mj_wk = mⱼ / ρⱼ * kernel_value( k, h_inv, xᵢ, xⱼ)

    map(Aⱼ) do Aj
        mj_wk * Aj
    end
end


"""
    𝒜(k::AbstractSPHKernel, h_inv, xᵢ, xⱼ, Aⱼ, mⱼ, ρⱼ )

Compute the contribution of particle `j` to the SPH quantity `A` for particle `i`.
Based on positions `xᵢ` and `xⱼ`.

``\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} W(\\vec{x}_i - \\vec{x}_j, h_i)``
"""
const 𝒜 = kernel_quantity
