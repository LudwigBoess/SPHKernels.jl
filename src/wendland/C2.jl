struct WendlandC2{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

struct WendlandC2_1D{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

"""
    WendlandC2(T::DataType=Float64, n_neighbours::Integer=216)

Set up a `WendlandC2` kernel for a given DataType `T`.
"""
function WendlandC2(T::DataType=Float64, dim::Integer=3)

    if dim == 1
        return WendlandC2_1D{T}(1, 5/4)
    elseif dim == 2 
        norm = T(7/π)
    elseif dim == 3
        norm = T(21/2π)
    else
        error("WendlandC2 not defined for $dim dimensions!")
    end

    return WendlandC2{T}(dim, norm)

end

"""
    WendlandC2(dim::Integer)

Define `WendlandC2` kernel with dimension `dim` for the native `DataType` of the OS.
"""
WendlandC2(dim::Integer) = WendlandC2(typeof(1.0), dim)

"""
    kernel_value(kernel::WendlandC2_1D{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::WendlandC2_1D{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv
        t1 = 1 - u
        t3 = t1*t1*t1
        return ( t3 * ( 1 + 3u )) * n |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::WendlandC2_1D{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^2
        t1 = 1 - u
        return ( -12u * t1  ) * n |> T
    else
        return 0 |> T
    end

end


"""
    kernel_value(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^kernel.dim
        t1 = 1 - u
        t4 = t1 * t1 * t1 * t1
        return ( t4 * ( 1 + 4u ) ) * n |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^kernel.dim * h_inv
        t1 = 1 - u
        t3 = t1 * t1 * t1
        return ( -20u * t3 ) * n |> T
    else
        return 0 |> T
    end

end


""" 
    bias_correction( kernel::Union{WendlandC2_1D{T}, WendlandC2{T}}, 
                     density::Real, m::Real, h_inv::Real,
                     n_neighbours::Integer ) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction( kernel::Union{WendlandC2_1D{T}, WendlandC2{T}}, 
                          density::Real, m::Real, h_inv::Real,
                          n_neighbours::Integer ) where T

    n = kernel.norm * h_inv^kernel.dim
    wc_correction = T(0.0294) * ( n_neighbours * T(0.01) )^T(-0.977) * m * n
    return density - wc_correction |> T
end