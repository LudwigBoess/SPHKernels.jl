struct Quintic{T} <: SPHKernel
    n_neighbours::Int64
    norm_1D::T
    norm_2D::T
    norm_3D::T 
end

"""
    Quintic(T::DataType=Float64, n_neighbours::Integer=64)

Set up a `Quintic` kernel for a given DataType `T`.
"""
Quintic(T::DataType=Float64, n_neighbours::Integer=216) = Quintic{T}(n_neighbours, 243.0/40.0, 15309.0/(478.0*π), 2187.0/(40.0*pi))
    

"""
    kernel_value_1D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

Evaluate quintic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_1D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_1D * h_inv

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
        return 0.0 |> T
    end

end

"""
    kernel_deriv_1D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the Quintic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_1D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_1D * h_inv^2

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
        return 0.0 |> T
    end

end

""" 
    bias_correction_1D(kernel::Quintic{T}, density::Real, m::Real, h_inv::Real) where T

Does not do anything for the BSplines. Implemented for stability.
"""
function bias_correction_1D(kernel::Quintic{T}, density::Real, m::Real, h_inv::Real) where T
    return density
end

"""
    kernel_value_2D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

Evaluate quintic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_2D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_2D * h_inv^2

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
        return 0.0 |> T
    end

end

"""
    kernel_deriv_2D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the Quintic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_2D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_2D * h_inv^3

    if u < 1/3
        u_m1  = 1 - u 
        u_m23 = T(2/3) - u
        u_m13 = T(1/3) - u
        return ( -5u_m1*u_m1*u_m1*u_m1 + 
                 30u_m23*u_m23*u_m23*u_m23 - 
                 75u_m13*u_m13*u_m13*u_m13 ) * n |> T
    elseif u < 2/3
        u_m1  = 1 - u 
        u_m23 = T(2/3) - u 
        return ( -5u_m1*u_m1*u_m1*u_m1 + 
                 30u_m23*u_m23*u_m23*u_m23 - 75 ) * n |> T
    elseif u < 1
        u_m1  = 1 - u
        return -5u_m1*u_m1*u_m1*u_m1 * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_2D(kernel::Quintic{T}, density::Real, m::Real, h_inv::Real) where T

Does not do anything for the BSplines. Implemented for stability.
"""
function bias_correction_2D(kernel::Quintic{T}, density::Real, m::Real, h_inv::Real) where T
    return density
end

"""
    kernel_value_3D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

Evaluate quintic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_3D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_3D * h_inv^3

    if u < 1/3
        u_m1  = 1 - u 
        u_m23 = T(2/3) - u 
        u_m13 = T(1/3) - u 
        return (  u_m1*u_m1*u_m1*u_m1*u_m1 - 
                 6u_m23*u_m23*u_m23*u_m23*u_m23 + 
                15u_m13*u_m13*u_m13*u_m13*u_m13 ) * n |> T
    elseif u < 2/3
        u_m1  = 1 - u 
        u_m23 = T(2/3) - u 
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 
                6u_m23*u_m23*u_m23*u_m23*u_m23 ) * n |> T
    elseif u < 1
        u_m1 = 1 - u 
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 ) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_3D(kernel::Quintic, u::Real, h_inv::Real)

Evaluate the derivative of the Quintic spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_3D(kernel::Quintic{T}, u::Real, h_inv::Real) where T

    n = kernel.norm_3D * h_inv^4

    if u < 1/3
        u_m1  = 1 - u 
        u_m23 = T(2/3) - u
        u_m13 = T(1/3) - u
        return ( -5u_m1*u_m1*u_m1*u_m1 + 
                30u_m23*u_m23*u_m23*u_m23 - 
                75u_m13*u_m13*u_m13*u_m13 ) * n |> T
    elseif u < 2/3
        u_m1  = 1 - u 
        u_m23 = T(2/3) - u 
        return ( -5u_m1*u_m1*u_m1*u_m1 + 
                 30u_m23*u_m23*u_m23*u_m23 - 75 ) * n |> T
    elseif u < 1
        u_m1  = 1 - u 
        return ( -5u_m1*u_m1*u_m1*u_m1 ) * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_3D(kernel::Quintic{T}, density::Real, m::Real, h_inv::Real) where T

Does not do anything for the BSplines. Implemented for stability.
"""
function bias_correction_3D(kernel::Quintic{T}, density::Real, m::Real, h_inv::Real) where T
    return density |> T
end