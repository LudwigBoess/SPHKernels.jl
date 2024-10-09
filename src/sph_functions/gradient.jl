
"""
    kernel_gradient( k::AbstractSPHKernel, h_inv::Real, xᵢ::T, xⱼ::T ) where T

Computes the gradient of the kernel `k` at the position of the neighbour `xⱼ`. 

``∇W(x_{ij}, h_i) = \\frac{dW}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}`` 
"""
function kernel_gradient( k::AbstractSPHKernel, h_inv::Real, xᵢ::T, xⱼ::T ) where T
    r  = get_r(xᵢ, xⱼ)

    dwk_r = d𝒲(k, r*h_inv, h_inv) / r

    map(xᵢ, xⱼ) do i, j
        dwk_r * (i - j)
    end
end

"""
    kernel_gradient( k::AbstractSPHKernel, r::T1, h_inv::T1, Δx::T2) where {T1,T2}

Computes the gradient of the kernel `k` at the distance `r` along the distance vector `Δx` of the neighbour `j`. 

``∇W(x_{ij}, h_i) = \\frac{dW}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}`` 
"""
function kernel_gradient( k::AbstractSPHKernel, r::T1, h_inv::T1, Δx::T2) where {T1,T2}

    dwk_r = d𝒲(k, r*h_inv, h_inv) / r

    map(Δx) do dx
        dwk_r * dx
    end

end


"""
    ∇𝒲( k::AbstractSPHKernel, h_inv::T1, xᵢ::T2, xⱼ::T2 ) where {T1<:Real,T2}

Computes the gradient of the kernel `k` at the position of the neighbour `xⱼ`. 
Compact notation of [`kernel_gradient`](@ref).

``∇W(x_{ij}, h_i) = \\frac{dW}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}`` 
"""
∇𝒲( k::AbstractSPHKernel, h_inv::T1, xᵢ::T2, xⱼ::T2 ) where {T1<:Real,T2} = kernel_gradient(k, h_inv, xᵢ, xⱼ)

"""
    ∇𝒲( k::AbstractSPHKernel, h_inv::Real, xᵢ::Union{Real, Vector{<:Real}}, xⱼ::Union{Real, Vector{<:Real}} )

Computes the gradient of the kernel `k` at the position of the neighbour `j`. 
Based on Euclidean distance `r` and distance vector `Δx` between the particles. 
Useful if many quantities need to be computed for the same particle pair.
Compact notation of [`kernel_gradient`](@ref).

``∇W(x_{ij}, h_i) = \\frac{dW}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}`` 
"""
∇𝒲( k::AbstractSPHKernel, r::T1, h_inv::T1, Δx::T2) where {T1<:Real,T2} = kernel_gradient(k, r, h_inv, Δx )


"""
    quantity_gradient(k::AbstractSPHKernel, h_inv::T1, 
                      xᵢ::T2, xⱼ::T2, Aⱼ::T2,
                      mⱼ::T1, ρⱼ::T1 ) where {T1<:Real, T2}

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.
Based on positions `xᵢ` and `xⱼ`.

``∇\\vec{A}_i(x) ≈ \\sum_j \\frac{m_j}{\\rho_j} \\vec{A}_j \\: ∇W(||\\vec{x}_i - \\vec{x}_j||, h_i)``
"""
function quantity_gradient( k::AbstractSPHKernel, h_inv::T1, 
                            xᵢ::T2, xⱼ::T2, Aⱼ::T2,
                            mⱼ::T1, ρⱼ::T1 ) where {T1<:Real, T2}
                  
    r = get_r(xᵢ, xⱼ)

    mj_dwk_r = mⱼ / (ρⱼ * r) * d𝒲(k, r * h_inv, h_inv)

    map(xᵢ, xⱼ, Aⱼ) do xi, xj, Aj
        mj_dwk_r * (xi - xj) * Aj
    end

end


"""
    quantity_gradient(k::AbstractSPHKernel,
                      r::T1, h_inv::T1,
                      Δx::T2, Aⱼ::T2,
                      mⱼ::T1, ρⱼ::T1) where {T1<:Real,T2}

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.
Based on Euclidean distance `r` and distance vector `Δx` between the particles. 
Useful if many quantities need to be computed for the same particle pair.

``∇\\vec{A}_i(x) ≈ \\sum_j \\frac{m_j}{\\rho_j} \\vec{A}_j \\: ∇W(||\\vec{x}_i - \\vec{x}_j||, h_i)``
"""
function quantity_gradient(k::AbstractSPHKernel,
                           r::T1, h_inv::T1,
                           Δx::T2, Aⱼ::T2,
                           mⱼ::T1, ρⱼ::T1) where {T1<:Real,T2}

    mj_dwk_r = mⱼ / (ρⱼ * r) * d𝒲(k, r * h_inv, h_inv)

    map(Δx, Aⱼ) do dx, Aj
        mj_dwk_r * dx * Aj
    end

end


"""
    ∇𝒜( k::AbstractSPHKernel, h_inv, xᵢ, xⱼ, Aⱼ, mⱼ, ρⱼ)

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.
Compact notation of [`quantity_gradient`](@ref).

``∇\\vec{A}_i(x) ≈ \\sum_j \\frac{m_j}{\\rho_j} \\vec{A}_j \\: ∇W(||\\vec{x}_i - \\vec{x}_j||, h_i)``
"""
∇𝒜( k::AbstractSPHKernel, h_inv, xᵢ, xⱼ, Aⱼ, mⱼ, ρⱼ) = quantity_gradient( k, h_inv, xᵢ, xⱼ, Aⱼ, mⱼ, ρⱼ)
