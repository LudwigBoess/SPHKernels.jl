struct WendlandC6{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

struct WendlandC6_1D{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

"""
    WendlandC6(T::DataType=Float64, dim::Integer=3)

Set up a `WendlandC6` kernel for a given DataType `T`.
"""
function WendlandC6(T::DataType = Float64, dim::Integer = 3)

    if dim == 1
        return WendlandC6_1D{T}(1, 55 / 32)
    elseif dim == 2
        norm = T(78 / (7π))
    elseif dim == 3
        norm = T(1365 / (64π))
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

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_value(kernel::WendlandC6_1D{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        t6 = t1 * t1 * t1 * t1 * t1 * t1
        u2 = u * u
        return (t6 * t1 * (1 + 7u + 19u2 + 21u2 * u)) |> T
    else
        return 0 |> T
    end
end

"""
    kernel_value(kernel::WendlandC6_1D{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
kernel_value(kernel::WendlandC6_1D{T}, 
    u::Real, h_inv::Real) where {T} = T(kernel.norm * h_inv) * kernel_value(kernel, u)

"""
    kernel_deriv(kernel::WendlandC6_1D{T}, u::Real) where T

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_deriv(kernel::WendlandC6_1D{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        t6 = t1 * t1 * t1 * t1 * t1 * t1
        u2 = u * u
        return (-6t6 * u * (35u2 + 18u + 1)) |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv(kernel::WendlandC6_1D{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
kernel_deriv(kernel::WendlandC6_1D{T}, 
    u::Real, h_inv::Real) where {T} = T(kernel.norm * h_inv^2) *kernel_deriv(kernel, u)


"""
    kernel_value(kernel::WendlandC6{T}, u::Real) where T

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``, without normalisation.
"""
function kernel_value(kernel::WendlandC6{T}, u::Real) where {T}

    if u < 1
        t1 = 1 - u
        t1 = t1 * t1  # (1.0 - u)^2
        t1 = t1 * t1  # (1.0 - u)^4
        t1 = t1 * t1  # (1.0 - u)^8
        u2 = u * u
        return (t1 * (1 + 8u + 25u2 + 32u2 * u)) |> T
    else
        return 0 |> T
    end

end

"""
    kernel_value(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
kernel_value(kernel::WendlandC6{T}, 
    u::Real, h_inv::Real) where {T} = T(kernel.norm * h_inv^kernel.dim) * kernel_value(kernel, u)

"""
    kernel_deriv(kernel::WendlandC6, u::Real)

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}`` without normalisation.
"""
function kernel_deriv(kernel::WendlandC6{T}, u::Real) where {T}


    if u < 1
        t1 = 1 - u
        t7 = t1 * t1 * t1 * t1 * t1 * t1 * t1
        return (-22t7 * u * (16 * u^2 + 7u + 1)) |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv(kernel::WendlandC6, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
kernel_deriv(kernel::WendlandC6{T}, 
    u::Real, h_inv::Real) where {T} = T(kernel.norm * h_inv^kernel.dim * h_inv) * kernel_deriv(kernel, u)


""" 
    bias_correction( kernel::Union{WendlandC6_1D{T}, WendlandC6{T}}, 
                     density::Real, m::Real, h_inv::Real, 
                     n_neighbours::Integer ) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction(kernel::Union{WendlandC6_1D{T},WendlandC6{T}},
    density::Real, m::Real, h_inv::Real,
    n_neighbours::Integer) where {T}

    n = kernel.norm * h_inv^kernel.dim
    wc_correction = T(0.0116) * (n_neighbours * T(0.01))^T(-2.236) * m * n

    if wc_correction < T(0.2) * density
        density -= wc_correction
    end

    return density |> T
end