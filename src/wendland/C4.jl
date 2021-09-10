struct WendlandC4{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

struct WendlandC4_1D{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

"""
    WendlandC4(T::DataType=Float64, n_neighbours::Integer=216)

Set up a `WendlandC4` kernel for a given DataType `T`.
"""
function WendlandC4(T::DataType=Float64, dim::Integer=3)

    if dim == 1
        return WendlandC4_1D{T}(1, 3/2)
    elseif dim == 2 
        norm = T(9/π)
    elseif dim == 3
        norm = T(495/32π)
    else
        error("WendlandC4 not defined for $dim dimensions!")
    end

    return WendlandC4{T}(dim, norm)

end

"""
    WendlandC4(dim::Integer)

Define `WendlandC4` kernel with dimension `dim` for the native `DataType` of the OS.
"""
WendlandC4(dim::Integer) = WendlandC4(typeof(1.0), dim)

"""
    kernel_value(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::WendlandC4_1D{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv
        t1 = 1 - u
        t5 = t1*t1*t1*t1*t1
        return ( t5 * (1 + 5u + 8u^2 )) * n |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv(kernel::WendlandC4_1D{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::WendlandC4_1D{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^2
        t1 = 1 - u
        t4 = t1*t1*t1*t1
        return ( -14u * t4 - 56u^2 * t4 ) * n |> T
    else
        return 0 |> T
    end

end


"""
    (kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^kernel.dim
        t1 = 1 - u
        t5 = t1*t1*t1*t1*t1
        return ( t1 * t5 * ( 1 + 6u + 35/3 * u^2 ) ) * n |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^kernel.dim * h_inv
        t1 = 1 - u
        t5 = t1*t1*t1*t1*t1
        return ( -288/3 * t5 * u^2 - 56/3 * u * t5 ) * n |> T
    else
        return 0 |> T
    end

end


""" 
    bias_correction( kernel::Union{WendlandC4_1D{T}, WendlandC4{T}}, 
                     density::Real, m::Real, h_inv::Real, 
                     n_neighbours::Int64 ) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction( kernel::Union{WendlandC4_1D{T}, WendlandC4{T}}, 
                          density::Real, m::Real, h_inv::Real, 
                          n_neighbours::Integer ) where T

    n = kernel.norm * h_inv^kernel.dim
    wc_correction = T(0.01342) * ( n_neighbours * T(0.01) )^T(-1.579) * m * n
    
    return density - wc_correction |> T
end