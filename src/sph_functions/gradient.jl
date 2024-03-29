
"""
    kernel_gradient( k::AbstractSPHKernel, h_inv::Real, 
                     xᵢ::Union{Real, Vector{<:Real}}, 
                     xⱼ::Union{Real, Vector{<:Real}} )

Computes the gradient of the kernel `k` at the position of the neighbour `xⱼ`. 

``∇W(x_{ij}, h_i) = \\frac{dW}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}`` 
"""
function kernel_gradient( k::AbstractSPHKernel, h_inv::Real, 
                          xᵢ::Vector{<:Real}, 
                          xⱼ::Vector{<:Real} )
    r  = get_r(xᵢ, xⱼ)

    dwk_r = d𝒲(k, r*h_inv, h_inv) / r

    map(xᵢ, xⱼ) do i, j
        dwk_r * (i - j)
    end
end

"""
    kernel_gradient( k::AbstractSPHKernel, r::Real, h_inv::Real, 
                     Δx::Vector{<:Real})

Computes the gradient of the kernel `k` at the distance `r` along the distance vector `Δx` of the neighbour `j`. 

``∇W(x_{ij}, h_i) = \\frac{dW}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}`` 
"""
function kernel_gradient( k::AbstractSPHKernel, r::Real, h_inv::Real, 
                          Δx::Union{Real, Vector{<:Real}})

    dwk_r = d𝒲(k, r*h_inv, h_inv) / r

    map(Δx) do dx
        dwk_r * dx
    end

end


"""
    kernel_gradient( k::AbstractSPHKernel, h_inv::Real, 
                     xᵢ::Union{Real, Vector{<:Real}}, 
                     xⱼ::Union{Real, Vector{<:Real}} )

Computes the gradient of the kernel `k` at the position of the neighbour `xⱼ`. 

``∇W(x_{ij}, h_i) = \\frac{dW}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}`` 
"""
function kernel_gradient( k::AbstractSPHKernel, h_inv::Real, 
                          xᵢ::Real,  xⱼ::Real, )
    
    r  = get_r(xᵢ, xⱼ)

    dwk_r = d𝒲(k, r*h_inv, h_inv) / r
    
    map(xᵢ, xⱼ) do i, j
        dwk_r * (i - j)
    end
end

"""
    ∇𝒲( k::AbstractSPHKernel, h_inv::Real, xᵢ::Union{Real, Vector{<:Real}}, xⱼ::Union{Real, Vector{<:Real}} )

Computes the gradient of the kernel `k` at the position of the neighbour `xⱼ`. 
Compact notation of [`kernel_gradient`](@ref).

``∇W(x_{ij}, h_i) = \\frac{dW}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}`` 
"""
∇𝒲( k::AbstractSPHKernel, h_inv::Real, 
     xᵢ::Union{Real, Vector{<:Real}}, 
     xⱼ::Union{Real, Vector{<:Real}} ) = kernel_gradient(k, h_inv, xᵢ, xⱼ)

"""
    ∇𝒲( k::AbstractSPHKernel, h_inv::Real, xᵢ::Union{Real, Vector{<:Real}}, xⱼ::Union{Real, Vector{<:Real}} )

Computes the gradient of the kernel `k` at the position of the neighbour `j`. 
Based on Euclidean distance `r` and distance vector `Δx` between the particles. 
Useful if many quantities need to be computed for the same particle pair.
Compact notation of [`kernel_gradient`](@ref).

``∇W(x_{ij}, h_i) = \\frac{dW}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}`` 
"""
∇𝒲( k::AbstractSPHKernel, r::Real, h_inv::Real, 
     Δx::Union{Real, Vector{<:Real}}) = kernel_gradient(k, r, h_inv, Δx )


"""
    quantity_gradient( k::AbstractSPHKernel, h_inv::Real, 
                       xᵢ::Union{Real, Vector{<:Real}},   
                       xⱼ::Union{Real, Vector{<:Real}},
                       Aⱼ::Union{Real, Vector{<:Real}},
                       mⱼ::Real,             ρⱼ::Real ) 

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.
Based on positions `xᵢ` and `xⱼ`.

``∇\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} ∇W(||\\vec{x}_i - \\vec{x}_j||, h_i)``
"""
function quantity_gradient( k::AbstractSPHKernel, h_inv::Real, 
                            xᵢ::Union{Real, Vector{<:Real}},   
                            xⱼ::Union{Real, Vector{<:Real}},
                            Aⱼ::Union{Real, Vector{<:Real}},
                            mⱼ::Real,             ρⱼ::Real ) 
                  
    r = get_r(xᵢ, xⱼ)

    mj_dwk_r = mⱼ / (ρⱼ * r) * d𝒲(k, r * h_inv, h_inv)

    map(xᵢ, xⱼ, Aⱼ) do xi, xj, Aj
        mj_dwk_r * (xi - xj) * Aj
    end

end


"""
    ∇𝒜( k::AbstractSPHKernel, h_inv::Real, 
        xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
        Aⱼ::Vector{<:Real},   
        mⱼ::Real=1,           ρⱼ::Real=1 )

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.
Compact notation of [`quantity_gradient`](@ref).

``∇\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} ∇W(||\\vec{x}_i - \\vec{x}_j||, h_i)``
"""
∇𝒜( k::AbstractSPHKernel, h_inv::Real, 
     xᵢ::Union{Real, Vector{<:Real}},   xⱼ::Union{Real, Vector{<:Real}},
     Aⱼ::Union{Real, Vector{<:Real}},
     mⱼ::Real,             ρⱼ::Real  ) = quantity_gradient( k, h_inv, 
                                                            xᵢ, xⱼ, 
                                                            Aⱼ, mⱼ, ρⱼ)

"""
    quantity_gradient( k::AbstractSPHKernel, 
                       r::Real,  h_inv::Real, 
                       Δx::Union{Real, Vector{<:Real}},
                       Aⱼ::Union{Real, Vector{<:Real}},
                       mⱼ::Real, ρⱼ::Real ) 

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.
Based on Euclidean distance `r` and distance vector `Δx` between the particles. 
Useful if many quantities need to be computed for the same particle pair.

``∇\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} ∇W(||\\vec{x}_i - \\vec{x}_j||, h_i)``
"""
function quantity_gradient( k::AbstractSPHKernel, 
                            r::Real,  h_inv::Real, 
                            Δx::Union{Real, Vector{<:Real}},
                            Aⱼ::Union{Real, Vector{<:Real}},
                            mⱼ::Real, ρⱼ::Real ) 

    mj_dwk_r = mⱼ / (ρⱼ * r) * d𝒲(k, r * h_inv, h_inv)

    map(Δx, Aⱼ) do dx, Aj
        mj_dwk_r * dx * Aj
    end

end


"""
    ∇𝒜( k::AbstractSPHKernel, 
         r::Real,  h_inv::Real, 
         Δx::Union{Real, Vector{<:Real}},
         Aⱼ::Union{Real, Vector{<:Real}},
         mⱼ::Real, ρⱼ::Real )

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.
Based on Euclidean distance `r` and distance vector `Δx` between the particles. 
Useful if many quantities need to be computed for the same particle pair.

``∇\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} ∇W(||\\vec{x}_i - \\vec{x}_j||, h_i)``
"""
∇𝒜( k::AbstractSPHKernel, 
     r::Real,  h_inv::Real, 
     Δx::Union{Real, Vector{<:Real}},
     Aⱼ::Vector{<:Real},
     mⱼ::Real, ρⱼ::Real )  = quantity_gradient( k, r, h_inv, Δx, Aⱼ, mⱼ, ρⱼ)


"""
    ∇𝒜( k::AbstractSPHKernel, 
         r::Real,  h_inv::Real, 
         Δx::Union{Real, Vector{<:Real}},
         Aⱼ::Union{Real, Vector{<:Real}},
         mⱼ::Real, ρⱼ::Real )

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.
Based on Euclidean distance `r` and distance vector `Δx` between the particles. 
Useful if many quantities need to be computed for the same particle pair.

``∇\\vec{A}_i(x) ≈ \\sum_j m_j \\frac{\\vec{A}_j}{\\rho_j} ∇W(||\\vec{x}_i - \\vec{x}_j||, h_i)``
"""
∇𝒜( k::AbstractSPHKernel, 
     r::Real,  h_inv::Real, 
     Δx::Union{Real, Vector{<:Real}},
     Aⱼ::Real,
     mⱼ::Real, ρⱼ::Real )  = Aⱼ * mⱼ / (ρⱼ * r) * d𝒲(k, r * h_inv, h_inv)
