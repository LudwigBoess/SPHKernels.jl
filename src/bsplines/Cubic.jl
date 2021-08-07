struct Cubic{T} <: SPHKernel
    n_neighbours::Int64
    norm_1D::T
    norm_2D::T
    norm_3D::T
end

"""
    Cubic(T::DataType=Float64, n_neighbours::Integer=64)

Set up a `Cubic` kernel for a given DataType `T`.
"""
Cubic(T::DataType=Float64, n_neighbours::Integer=64) = Cubic{T}(n_neighbours, 4.0/3.0, 40.0/7π, 8.0/π)

"""
    kernel_value_1D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

Evaluate cubic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_1D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_1D * h_inv

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
    kernel_deriv_1D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the Cubic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_1D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_1D * h_inv^2

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
    bias_correction_1D(kernel::Cubic{T}, density::Real, m::Real, h_inv::Real) where T

Does not do anything for the BSplines. Implemented for stability.
"""
function bias_correction_1D(kernel::Cubic{T}, density::Real, m::Real, h_inv::Real) where T
    return density
end

"""
    kernel_value_2D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

Evaluate cubic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_2D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_2D * h_inv^2

    if u < 0.5
        return ( 1 + 6 * ( u - 1 ) * u^2) * n |> T
    elseif u < 1.0
        u_m1 = 1 - u
        return  2u_m1^3 * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_2D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the Cubic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_2D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_2D * h_inv^3

    if u < 0.5
        return ( u * (18u - 12 )) * n |> T
    elseif u < 1.0
        u_m1 = 1 - u 
        return  -6u_m1^2 * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_2D(kernel::Cubic, density::Real, m::Real, h_inv::Real)

Does not do anything for the BSplines. Implemented for stability.
"""
function bias_correction_2D(kernel::Cubic{T}, density::Real, m::Real, h_inv::Real) where T
    return density
end


"""
    kernel_value_3D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

Evaluate cubic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_3D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_3D * h_inv^3
    
    if u < 0.5
        return ( 1 + 6 * ( u - 1 ) * u^2) * n |> T
    elseif u < 1.0
        u_m1 = 1 - u 
        return 2u_m1^3 * n |> T
    else
        return 0.0 |> T
    end

end


"""
    kernel_deriv_3D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the Cubic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_3D(kernel::Cubic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_3D * h_inv^2 * h_inv^2

    if u < 0.5
        return ( u * (18u - 12 )) * n |> T
    elseif u < 1.0
        u_m1 =  1 - u
        return -6u_m1^2 * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_3D(kernel::Cubic{T}, density::Real, m::Real, h_inv::Real) where T

Does not do anything for the BSplines. Implemented for stability.
"""
function bias_correction_3D(kernel::Cubic{T}, density::Real, m::Real, h_inv::Real)  where T
    return density |> T
end