struct WendlandC4{T} <: WendlandKernel
    dim::Int8
    norm::T
end

struct WendlandC4_1D{T} <: WendlandKernel
    dim::Int8
    norm::T
end

"""
    WendlandC4(T::DataType=Float64, dim::Integer=3)

Set up a `WendlandC4` kernel for a given DataType `T`.
"""
function WendlandC4(T::DataType = Float64, dim::Integer = 3)

    if dim == 1
        return WendlandC4_1D{T}(1, 3 / 2)
    elseif dim == 2
        norm = T(9 / π)
    elseif dim == 3
        norm = T(495 / 32π)
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
    kernel_value(kernel::WendlandC4{T}, u::Real) where T

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_value(kernel::WendlandC4_1D{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        t4 = t1 * t1    # (1-u)^2
        t4 *= t4        # (1-u)^4
        return t4 * t1 * (1 + 5u + 8u^2)
    else
        return zero(T)
    end

end


"""
    kernel_deriv(kernel::WendlandC4_1D{T}, u::Real) where T

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_deriv(kernel::WendlandC4_1D{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        t1 *= t1    # (1-u)^2
        t1 *= t1    # (1-u)^4
        return -14u * t1 - 56u^2 * t1 |> T
    else
        return zero(T)
    end

end


"""
    kernel_value(kernel::WendlandC4{T}, u::Real) where T

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_value(kernel::WendlandC4{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        t1 *= t1       # (1-u)^2
        t4 = t1 * t1   # (1-u)^4
        return t1 * t4 * (1 + 6u + 35 / 3 * u*u) |> T
    else
        return zero(T)
    end

end


"""
    kernel_deriv(kernel::WendlandC4{T}, u::Real) where T

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_deriv(kernel::WendlandC4{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        t5 = t1 * t1   # (1-u)^2
        t5 *= t5       # (1-u)^4
        t5 *= t1       # (1-u)^5
        return (-288 / 3 * t5 * u^2 - 56 / 3 * u * t5) |> T
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
function bias_correction(kernel::Union{WendlandC4_1D{T},WendlandC4{T}},
                        density::Real, m::Real, h_inv::Real,
                        n_neighbours::Integer) where {T}

    n = kernel_norm(kernel, h_inv)
    wc_correction = T(0.01342) * (n_neighbours * T(0.01))^T(-1.579) * m * n

    return density - wc_correction |> T
end