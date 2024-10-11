
"""
    density_estimate(k::AbstractSPHKernel, m::Real, r::Real, h_inv::Real)

Contribution of particle j to the density estimate of particle i.
"""
function density_estimate(k::AbstractSPHKernel, m::T, r::T, h_inv::T) where T
    return m * 𝒲(k, r * h_inv, h_inv)
end

"""
    ρⱼ(k::AbstractSPHKernel, m::Real, r::Real, h_inv::Real)

Contribution of particle j to the density estimate of particle i.
See also [`density_estimate`](@ref).
"""
ρⱼ(k::AbstractSPHKernel, m::T, r::T, h_inv::T) where T = density_estimate(k, m, r, h_inv)
