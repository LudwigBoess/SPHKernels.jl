struct Quintic{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

"""
    Quintic(T::DataType=Float64, n_neighbours::Integer=64)

Set up a `Quintic` kernel for a given DataType `T`.
"""    
function Quintic(T::DataType=Float64, dim::Integer=3)

    if dim == 1
        norm = 243/40
    elseif dim == 2 
        norm = 15309/478π
    elseif dim == 3
        norm = 2187/40π
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
    kernel_value(kernel::Quintic{T}, u::Real, h_inv::Real) where T

Evaluate quintic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::Quintic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm * h_inv^kernel.dim

    if u < 1/3
        u_m1  = 1 - u 
        u_m23 = T(2/3) - u 
        u_m13 = T(1/3) - u 
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 
                6u_m23*u_m23*u_m23*u_m23*u_m23 + 
                15u_m13*u_m13*u_m13*u_m13*u_m13 ) * n |> T
    elseif u < 2/3
        u_m1  = 1 - u
        u_m23 = T(2/3) - u 
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 
                6u_m23*u_m23*u_m23*u_m23*u_m23 ) * n |> T
    elseif u < 1
        u_m1  = 1 - u 
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 ) * n |> T
    else
        return 0 |> T
    end

end

"""
    kernel_deriv(kernel::Quintic{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the Quintic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::Quintic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm * h_inv^(kernel.dim + 1)

    if u < 1/3
        u_m1  = 1 - u 
        u_m23 = T(2/3) - u 
        u_m13 = T(1/3) - u 
        return ( -5u_m1*u_m1*u_m1*u_m1 + 
                 30u_m23*u_m23*u_m23*u_m23  - 
                 75u_m13*u_m13*u_m13*u_m13 ) * n |> T
    elseif u < 2/3
        u_m1  = 1 - u 
        u_m23 = T(2/3) - u
        return ( -5u_m1*u_m1*u_m1*u_m1 + 
                 30u_m23*u_m23*u_m23*u_m23 - 75 ) * n |> T
    elseif u < 1
        u_m1 =  1 - u 
        return -5u_m1*u_m1*u_m1*u_m1 * n |> T
    else
        return 0 |> T
    end

end

""" 
    bias_correction( kernel::Quintic{T}, 
                     density::Real, m::Real, h_inv::Real, 
                     n_neighbours::Integer ) where T

Does not do anything for the BSplines. Implemented for stability.
"""
function bias_correction( kernel::Quintic{T}, 
                          density::Real, m::Real, h_inv::Real, 
                          n_neighbours::Integer ) where T
    return density
end
