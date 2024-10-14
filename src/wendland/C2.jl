struct WendlandC2{T} <: WendlandKernel
    dim::Int8
    norm::T
end

struct WendlandC2_1D{T} <: WendlandKernel
    dim::Int8
    norm::T
end

"""
    WendlandC2(T::DataType=Float64, dim::Integer=3)

Set up a `WendlandC2` kernel for a given DataType `T`.
"""
function WendlandC2(T::DataType = Float64, dim::Integer = 3)

    if dim == 1
        return WendlandC2_1D{T}(1, 5 / 4)
    elseif dim == 2
        norm = T(7 / π)
    elseif dim == 3
        norm = T(21 / 2π)
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
    kernel_value(kernel::WendlandC2_1D{T}, u::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_value(kernel::WendlandC2_1D{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        return t1 * t1 * t1 * (1 + 3u)  |> T
    else
        return zero(T)
    end

end


"""
    kernel_deriv_1D(kernel::WendlandC2{T}, u::Real) where T

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_deriv(kernel::WendlandC2_1D{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        return -12u * t1 * t1 * t1 |> T
    else
        return zero(T)
    end

end


"""
    kernel_value(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_value(kernel::WendlandC2{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        t1 *= t1    # (1-u)^2
        t1 *= t1    # (1-u)^4
        return t1 * (1 + 4u) |> T
    else
        return zero(T)
    end

end


"""
    kernel_deriv(kernel::WendlandC2{T}, u::Real) where T

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_deriv(kernel::WendlandC2{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        return -20u * t1 * t1 * t1 |> T
    else
        return zero(T)
    end

end


""" 
    bias_correction( kernel::Union{WendlandC2_1D{T}, WendlandC2{T}}, 
                     density::Real, m::Real, h_inv::Real,
                     n_neighbours::Integer ) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction(kernel::Union{WendlandC2_1D{T},WendlandC2{T}},
    density::Real, m::Real, h_inv::Real,
    n_neighbours::Integer) where {T}

    n = kernel_norm(kernel, h_inv)
    wc_correction = T(0.0294) * (n_neighbours * T(0.01))^T(-0.977) * m * n
    return density - wc_correction |> T
end