
function kernel_rotation(k::AbstractSPHKernel, h_inv::Real, 
                         part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real},
                         part_Q::Vector{<:Real}, neighbour_Q::Vector{<:Real})

    Δx = part_pos[1] - neighbour_pos[1]
    Δy = part_pos[2] - neighbour_pos[2]
    Δz = part_pos[3] - neighbour_pos[3]

    ΔQx = part_Q[1] - neighbour_Q[1]
    ΔQy = part_Q[2] - neighbour_Q[2]
    ΔQz = part_Q[3] - neighbour_Q[3]

    r = √(Δx^2 + Δy^2 + Δz^2)
    u = r * h_inv

    dWk_r = d𝒲₃(k, u, h_inv) / r
    
    return  dWk_r * ( Δy * ΔQz - Δz * ΔQy), 
            dWk_r * ( Δz * ΔQx - Δx * ΔQz), 
            dWk_r * ( Δx * ΔQy - Δy * ΔQx)

end


∇xQ( k::AbstractSPHKernel, h_inv::Real, 
      part_pos::Vector{<:Real}, neighbour_pos::Vector{<:Real},
      part_Q::Vector{<:Real}, neighbour_Q::Vector{<:Real} ) = kernel_rotation( k, h_inv, 
                                                                               part_pos, neighbour_pos, 
                                                                               part_Q, neighbour_Q)