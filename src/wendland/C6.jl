struct WendlandC6{T} <: SPHKernel
    dim::Int64
    norm::T
end

struct WendlandC6_1D{T} <: SPHKernel
    dim::Int64
    norm::T
end

"""
    WendlandC6(T::DataType=Float64, dim::Integer=3)

Set up a `WendlandC6` kernel for a given DataType `T`.
"""
function WendlandC6(T::DataType=Float64, dim::Integer=3)

    if dim == 1
        return WendlandC6_1D{T}(1, 55/32)
    elseif dim == 2 
        norm = T(78/(7π))
    elseif dim == 3
        norm = T(1365/(64π))
    else
        error("WendlandC6 not defined for $dim dimensions!")
    end

    return WendlandC6{T}(dim, norm)

end

"""
    WendlandC6(dim::Integer)

Define `WendlandC6` kernel with dimension `dim` for the native `DataType` of the OS.
"""
WendlandC6(dim::Integer) = WendlandC6(typeof(1.0), dim)

"""
    kernel_value(kernel::WendlandC6_1D{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::WendlandC6_1D{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv
        t1 = 1 - u
        t6 = t1*t1*t1*t1*t1*t1
        u2 = u*u
        return ( t6 * t1 * (1 + 7u + 19u2 + 21u2 * u)) * n |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv(kernel::WendlandC6_1D{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::WendlandC6_1D{T}, u::Real, h_inv::Real) where T

     if u < 1
        n = kernel.norm * h_inv^2
        t1 = 1 - u
        t6 = t1*t1*t1*t1*t1*t1
        u2 = u*u
        return ( -6t6 * u * (35u2 + 18u + 3 )) * n |> T
    else
        return 0 |> T
    end

end


"""
    kernel_value_2D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^kernel.dim
        u_m1 = 1 - u
        u_m1 = u_m1 * u_m1  # (1.0 - u)^2
        u_m1 = u_m1 * u_m1  # (1.0 - u)^4
        u_m1 = u_m1 * u_m1  # (1.0 - u)^8
        u2 = u*u
        return ( u_m1 * ( 1 + 8u + 25u2 + 32u2*u )) * n |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC6, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T


    if u < 1
        n = kernel.norm * h_inv^kernel.dim * h_inv
        u_m1 = 1 - u
        u_m1_7 = u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1
        return ( -22u_m1_7 * u * ( 16*u^2 + 7u + 1 )) * n |> T
    else
        return 0 |> T 
    end

end


""" 
    bias_correction( kernel::Union{WendlandC6_1D{T}, WendlandC6{T}}, 
                     density::Real, m::Real, h_inv::Real, 
                     n_neighbours::Integer ) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction(kernel::Union{WendlandC6_1D{T}, WendlandC6{T}}, 
                         density::Real, m::Real, h_inv::Real, 
                         n_neighbours::Integer) where T

    n = kernel.norm * h_inv^kernel.dim
    wc_correction = T(0.0116) * ( n_neighbours * T(0.01) )^T(-2.236) * m * n

    if wc_correction < T(0.2)*density
        density -= wc_correction
    end
    
    return density |> T
end