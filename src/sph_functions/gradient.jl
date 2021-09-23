
"""
    kernel_gradient( k::SPHKernel, h_inv::Real, 
                     xᵢ::Union{Real, Vector{<:Real}}, 
                     xⱼ::Union{Real, Vector{<:Real}} )

Computes the gradient of the kernel `k` at the position of the neighbour `xⱼ`. 

``∇W(x_{ij}, h_i) = \frac{dW}{dx}\vert_{x_j} \frac{Δx_{ij}}{||x_{ij}||} \frac{1}{h_i}`` 
"""
function kernel_gradient( k::SPHKernel, h_inv::Real, 
                          xᵢ::Union{Real, Vector{<:Real}}, 
                          xⱼ::Union{Real, Vector{<:Real}} )
    
    r = √( sum( (xᵢ .- xⱼ).^2 ) )

    q    = r*h_inv
    dWdq = d𝒲(k, q, h_inv)

    grad = map(xᵢ, xⱼ) do i, j
        d = i - j
        dWdq*(d/r)
    end
    
    return grad
end


"""
    ∇𝒲( k::SPHKernel, h_inv::Real, xᵢ::Union{Real, Vector{<:Real}}, xⱼ::Union{Real, Vector{<:Real}} )

Computes the gradient of the kernel `k` at the position of the neighbour `xⱼ`. 
Multiple dispatch version of [`kernel_gradient`](@ref).

``∇W(x_{ij}, h_i) = \frac{dW}{dx}\vert_{x_j} \frac{Δx_{ij}}{||x_{ij}||} \frac{1}{h_i}`` 
"""
∇𝒲( k::SPHKernel, h_inv::Real, 
     xᵢ::Union{Real, Vector{<:Real}}, 
     xⱼ::Union{Real, Vector{<:Real}} ) = kernel_gradient(k, h_inv, xᵢ, xⱼ)

"""
    quantity_gradient( k::AbstractSPHKernel, h_inv::Real, 
                       xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                       Aᵢ::Vector{<:Real},   Aⱼ::Vector{<:Real},
                       mⱼ::Real,           ρⱼ::Real )

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.

``∇\vec{A}_i(x) ≈ \sum_j m_j \frac{\vec{A}_j}{\rho_j} ∇W(||\vec{x}_i - \vec{x}_j||, h_i)``
"""
quantity_gradient( k::AbstractSPHKernel, h_inv::Real, 
                   xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
                   Aⱼ::Vector{<:Real},
                   mⱼ::Real,             ρⱼ::Real ) = mⱼ / ρⱼ * Aⱼ * ∇𝒲( k, h_inv, xᵢ, xⱼ)


"""
    ∇𝒜( k::AbstractSPHKernel, h_inv::Real, 
        xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
        Aᵢ::Vector{<:Real},   
        mⱼ::Real=1,           ρⱼ::Real=1 )

Compute the contribution of particle `j` to the gradient of the SPH quantity `A` for particle `i`.
Multiple dispatch version of [`quantity_gradient`](@ref).

``∇\vec{A}_i(x) ≈ \sum_j m_j \frac{\vec{A}_j}{\rho_j} ∇W(||\vec{x}_i - \vec{x}_j||, h_i)``
"""
∇𝒜( k::AbstractSPHKernel, h_inv::Real, 
     xᵢ::Vector{<:Real},   xⱼ::Vector{<:Real},
     Aⱼ::Vector{<:Real},
     mⱼ::Real,             ρⱼ::Real ) = quantity_gradient( k, h_inv, xᵢ, xⱼ)