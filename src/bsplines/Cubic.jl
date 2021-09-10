struct Cubic{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

"""
    Cubic(T::DataType=Float64, n_neighbours::Integer=64)

Set up a `Cubic` kernel for a given DataType `T` and dimensin `dim`.
"""
function Cubic(T::DataType=Float64, dim::Integer=3)

    if dim == 1
        norm = 4/3
    elseif dim == 2 
        norm = 40/7π 
    elseif dim == 3
        norm = 8/π
    else
        error("Cubic not defined for $dim dimensions!")
    end

    return Cubic{T}(dim, norm)

end

"""
    Cubic(dim::Integer)

Define `Cubic` kernel with dimension `dim` for the native `DataType` of the OS.
"""
Cubic(dim::Integer) = Cubic(typeof(1.0), dim)

"""
    kernel_value_1D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

Evaluate cubic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::Cubic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm * h_inv^kernel.dim

    if u < 0.5
        return ( 1 + 6 * ( u - 1 ) * u^2) * n |> T
    elseif u < 1
        u_m1 = 1 - u 
        return  2u_m1^3 * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv(kernel::Cubic{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the Cubic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::Cubic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm * h_inv^(kernel.dim+1)

    if u < 0.5
        return ( u * (18u - 12 )) * n |> T
    elseif u < 1
        u_m1 = 1 - u
        return -6u_m1^2 * n |> T
    else
        return 0.0 |> T
    end

end


""" 
    bias_correction( kernel::Cubic{T}, 
                     density::Real, m::Real, h_inv::Real, 
                     n_neighbours::Integer )  where T

Does not do anything for the BSplines. Implemented for stability.
"""
function bias_correction( kernel::Cubic{T}, 
                          density::Real, m::Real, h_inv::Real, 
                          n_neighbours::Integer )  where T
    return density |> T
end