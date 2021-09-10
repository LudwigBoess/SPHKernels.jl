struct WendlandC8{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

"""
    WendlandC8(T::DataType=Float64, n_neighbours::Integer=395)

Set up a `WendlandC8` kernel for a given DataType `T`.
"""
function WendlandC8(T::DataType=Float64, dim::Integer=3)

    if dim == 2
        norm = T(8/3π)
    elseif dim == 3
        norm = T(357/64π)
    else
        error("WendlandC8 not defined for $dim dimensions!")
    end

    return WendlandC8{T}(dim, norm)

end

"""
    WendlandC8(dim::Integer)

Define `WendlandC8` kernel with dimension `dim` for the native `DataType` of the OS.
"""
WendlandC8(dim::Integer) = WendlandC8(typeof(1.0), dim)

"""
    kernel_value(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC8 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^kernel.dim
        t1 = 1 - u
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return ( t1 * t9 * (5 + 50u + 210u2 + 450u2 * u + 429u2 * u2)) * n |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC8 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^kernel.dim * h_inv
        t1 = 1.0 - u
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return -26t9 * u * (231u2 * u + 159u2 + 45u + 5 ) * n |> T
    else
        return 0 |> T
    end

end

""" 
    bias_correction(kernel::WendlandC8, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction( kernel::WendlandC8{T}, 
                          density::Real, m::Real, h_inv::Real, 
                          n_neighbours::Integer ) where T
    return density |> T
end
