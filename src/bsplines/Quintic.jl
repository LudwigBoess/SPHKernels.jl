struct Quintic{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

"""
    Quintic(T::DataType=Float64, dim::Integer=3)

Set up a `Quintic` kernel for a given DataType `T`.
"""
function Quintic(T::DataType = Float64, dim::Integer = 3)

    if dim == 1
        norm = 243 / 40
    elseif dim == 2
        norm = 15309 / 478π
    elseif dim == 3
        norm = 2187 / 40π
    else
        error("Quintic not defined for $dim dimensions!")
    end

    return Quintic{T}(dim, norm)

end

"""
    Quintic(dim::Integer)

Define `Quintic` kernel with dimension `dim` for the native `DataType` of the OS.
"""
Quintic(dim::Integer) = Quintic(typeof(1.0), dim)

"""
    kernel_value(kernel::Quintic{T}, u::Real) where T

Evaluate quintic spline at position ``u = \\frac{x}{h}``, without normalisation.
"""
function kernel_value(kernel::Quintic{T}, u::Real) where {T}

    if u < 1 / 3
        u_m1 = 1 - u
        u_m23 = T(2 / 3) - u
        u_m13 = T(1 / 3) - u
        return (u_m1 * u_m1 * u_m1 * u_m1 * u_m1 -
                6u_m23 * u_m23 * u_m23 * u_m23 * u_m23 +
                15u_m13 * u_m13 * u_m13 * u_m13 * u_m13) |> T
    elseif u < 2 / 3
        u_m1 = 1 - u
        u_m23 = T(2 / 3) - u
        return (u_m1 * u_m1 * u_m1 * u_m1 * u_m1 -
                6u_m23 * u_m23 * u_m23 * u_m23 * u_m23) |> T
    elseif u < 1
        u_m1 = 1 - u
        return (u_m1 * u_m1 * u_m1 * u_m1 * u_m1) |> T
    else
        return 0 |> T
    end
end

"""
    kernel_value(kernel::Quintic{T}, u::Real, h_inv::Real) where T

Evaluate quintic spline at position ``u = \\frac{x}{h}``.
"""
kernel_value(kernel::Quintic{T}, u::Real, h_inv::Real) where {T} = 
    T(kernel.norm * h_inv^kernel.dim) * kernel_value(kernel, u)

"""
    kernel_deriv(kernel::Quintic{T}, u::Real) where T

Evaluate the derivative of the Quintic spline at position ``u = \\frac{x}{h}``, without normalisation.
"""
function kernel_deriv(kernel::Quintic{T}, u::Real) where {T}

    if u < 1 / 3
        u_m1 = 1 - u
        u_m23 = T(2 / 3) - u
        u_m13 = T(1 / 3) - u
        return (-5u_m1 * u_m1 * u_m1 * u_m1 +
                30u_m23 * u_m23 * u_m23 * u_m23 -
                75u_m13 * u_m13 * u_m13 * u_m13) |> T
    elseif u < 2 / 3
        u_m1 = 1 - u
        u_m23 = T(2 / 3) - u
        return (-5u_m1 * u_m1 * u_m1 * u_m1 +
                30u_m23 * u_m23 * u_m23 * u_m23) |> T
    elseif u < 1
        u_m1 = 1 - u
        return -5u_m1 * u_m1 * u_m1 * u_m1 |> T
    else
        return 0 |> T
    end
end

"""
    kernel_deriv(kernel::Quintic{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the Quintic spline at position ``u = \\frac{x}{h}``.
"""
kernel_deriv(kernel::Quintic{T}, u::Real, h_inv::Real) where {T} = 
    T(kernel.norm * h_inv^kernel.dim * h_inv) * kernel_deriv(kernel, u)

""" 
    bias_correction( kernel::Quintic{T}, 
                     density::Real, m::Real, h_inv::Real, 
                     n_neighbours::Integer ) where T

Does not do anything for the BSplines. Implemented for stability.
"""
function bias_correction(kernel::Quintic{T},
                        density::Real, m::Real, h_inv::Real,
                        n_neighbours::Integer) where {T}
    return density
end
