"""
    kernel_grad_1D(k::AbstractSPHKernel, h_inv::Real, part_pos::Real, neighbour_pos::Real)

Computes the 1D gradient of the kernel `k` at the position of the neighbour `neighbour_pos`. 

``∇𝒲(x_i) = \\frac{d𝒲}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}``
"""
function kernel_grad_1D(k::AbstractSPHKernel, h_inv::Real, part_pos::Real, neighbour_pos::Real)

    x_diff = part_pos - neighbour_pos

    r = √( x_diff^2 )
    u = r * h_inv

    return d𝒲₁(k, u, h_inv) / r * x_diff
end

""" 
    ∇𝒲₁(k::AbstractSPHKernel, h_inv::Real, part_pos::Real, neighbour_pos::Real)

Multiple dispatch version of [`kernel_grad_1D`](@ref).
"""
∇𝒲₁(k::AbstractSPHKernel, h_inv::Real, part_pos::Real, neighbour_pos::Real) = kernel_grad_1D(k, h_inv, part_pos, neighbour_pos)

""" 
    ∇𝒲₁(k::AbstractSPHKernel, h_inv::Real, part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real})

Multiple dispatch version of [`kernel_grad_1D`](@ref).
"""
∇𝒲₁(k::AbstractSPHKernel, h_inv::Real, part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real}) = kernel_grad_1D(k, h_inv, part_pos[1], neighbour_pos[1])



"""
    kernel_grad_2D(k::AbstractSPHKernel, h_inv::Real, part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real})

Computes the 2D gradient of the kernel `k` at the position of the neighbour `neighbour_pos`. 

``∇𝒲(x_i) = \\frac{d𝒲}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}``
"""
function kernel_grad_2D(k::AbstractSPHKernel, h_inv::Real, part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real})

    x_diff = part_pos[1] - neighbour_pos[1]
    y_diff = part_pos[2] - neighbour_pos[2]

    r = √(x_diff^2 + y_diff^2 )
    u = r * h_inv

    tmp = d𝒲₂(k, u, h_inv) / r
    Wgx = tmp*x_diff
    Wgy = tmp*y_diff

    return Wgx, Wgy
end

""" 
    ∇𝒲₂(k::AbstractSPHKernel, h_inv::Real, part_pos::Real, neighbour_pos::Real)

Multiple dispatch version of [`kernel_grad_2D`](@ref).
"""
∇𝒲₂(k::AbstractSPHKernel, h_inv::Real, part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real}) = kernel_grad_2D(k, h_inv, part_pos, neighbour_pos)


"""
    kernel_grad_3D(k::AbstractSPHKernel, h_inv::Real, part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real})

Computes the 3D gradient of the kernel `k` at the position of the neighbour `neighbour_pos`. 

``∇𝒲(x_i) = \\frac{d𝒲}{dx}\\vert_{x_j} \\frac{Δx_{ij}}{||x_{ij}||} \\frac{1}{h_i}``
"""
function kernel_grad_3D(k::AbstractSPHKernel, h_inv::Real, part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real})

    x_diff = part_pos[1] - neighbour_pos[1]
    y_diff = part_pos[2] - neighbour_pos[2]
    z_diff = part_pos[3] - neighbour_pos[3]

    r = √(x_diff^2 + y_diff^2 + z_diff^2)
    u = r * h_inv

    tmp = d𝒲₃(k, u, h_inv) / r
    Wgx = tmp*x_diff
    Wgy = tmp*y_diff
    Wgz = tmp*z_diff

    return Wgx, Wgy, Wgz
end


""" 
    ∇𝒲₃(k::AbstractSPHKernel, h_inv::Real, part_pos::Real, neighbour_pos::Real)

Multiple dispatch version of [`kernel_grad_3D`](@ref).
"""
∇𝒲₃(k::AbstractSPHKernel, h_inv::Real, part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real}) = kernel_grad_3D(k, h_inv, part_pos, neighbour_pos)

