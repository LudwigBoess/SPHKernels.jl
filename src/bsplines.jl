"""
    struct Cubic <: SPHKernel
        n_neighbours::Int64
        norm_1D::Float64
        norm_2D::Float64
        norm_3D::Float64
        function Cubic(n_neighbours::Int64=64)
            new(n_neighbours, 8.0/π, 8.0/π)
        end
    end

Datatype for cubic sph spline.
"""
struct Cubic <: SPHKernel
    n_neighbours::Int64
    norm_1D::Float64
    norm_2D::Float64
    norm_3D::Float64
    function Cubic(n_neighbours::Integer=64)
        new(n_neighbours, 4.0/3.0, 40.0/7π, 8.0/π)
    end
end

"""
    kernel_value_1D(kernel::Cubic, u::Real, h_inv::Real)

Evaluate cubic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_1D(kernel::Cubic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_1D * h_inv

    @fastmath if u < 0.5
        return ( 1.0 + 6.0 * ( u - 1.0 ) * u^2) * n
    elseif u < 1.0
        u_m1 = (1.0 - u )
        return ( 2.0 * u_m1^3 ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_1D(kernel::Cubic, u::Real, h_inv::Real)

Evaluate the derivative of the Cubic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_1D(kernel::Cubic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_1D * h_inv^2

    @fastmath if u < 0.5
        return ( u * (18.0 * u - 12.0 )) * n
    elseif u < 1.0
        u_m1 = ( 1.0 - u )
        return ( -6.0 * u_m1^2 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_1D(kernel::Cubic, density::Real, m::Real, h_inv::Real)

Does not do anything for the BSplines. Implemented for stability.
"""
@inline function bias_correction_1D(kernel::Cubic, density::Real, m::Real, h_inv::Real)  
    return density
end

"""
    kernel_value_2D(kernel::Cubic, u::Real, h_inv::Real)

Evaluate cubic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_2D(kernel::Cubic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_2D * h_inv^2

    @fastmath if u < 0.5
        return ( 1.0 + 6.0 * ( u - 1.0 ) * u^2) * n
    elseif u < 1.0
        u_m1 = (1.0 - u )
        return ( 2.0 * u_m1^3 ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_2D(kernel::Cubic, u::Real, h_inv::Real)

Evaluate the derivative of the Cubic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_2D(kernel::Cubic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_2D * h_inv^3

    @fastmath if u < 0.5
        return ( u * (18.0 * u - 12.0 )) * n
    elseif u < 1.0
        u_m1 = ( 1.0 - u )
        return ( -6.0 * u_m1^2 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_2D(kernel::Cubic, density::Real, m::Real, h_inv::Real)

Does not do anything for the BSplines. Implemented for stability.
"""
@inline function bias_correction_2D(kernel::Cubic, density::Real, m::Real, h_inv::Real)  
    return density
end


"""
    kernel_value_3D(kernel::Cubic, u::Real, h_inv::Real)

Evaluate cubic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_3D(kernel::Cubic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_3D * h_inv^3
    
    @fastmath if u < 0.5
        return ( 1.0 + 6.0 * ( u - 1.0 ) * u^2) * n
    elseif u < 1.0
        u_m1 = (1.0 - u )
        return ( 2.0 * u_m1^3 ) * n
    else
        return 0.
    end

end


"""
    kernel_deriv_3D(kernel::Cubic, u::Real, h_inv::Real)

Evaluate the derivative of the Cubic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_3D(kernel::Cubic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_3D * h_inv^4

    @fastmath if u < 0.5
        return ( u * (18.0 * u - 12.0 )) * n
    elseif u < 1.0
        u_m1 = ( 1.0 - u )
        return ( -6.0 * u_m1^2 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_3D(kernel::Cubic, density::Real, m::Real, h_inv::Real)

Does not do anything for the BSplines. Implemented for stability.
"""
@inline function bias_correction_3D(kernel::Cubic, density::Real, m::Real, h_inv::Real)  
    return density
end


"""
    struct Quintic <: SPHKernel
        n_neighbours::Int64
        norm_1D::Float64
        norm_2D::Float64
        norm_3D::Float64
    end

Datatype for quintic sph spline.
"""
struct Quintic <: SPHKernel
    n_neighbours::Int64
    norm_1D::Float64
    norm_2D::Float64
    norm_3D::Float64
    function Quintic(n_neighbours::Integer=216)
        new(n_neighbours, 243.0/40.0, 15309.0/(478.0*π), 2187.0/(40.0*pi))
    end
end

"""
    kernel_value_1D(kernel::Quintic, u::Real, h_inv::Real)

Evaluate quintic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_1D(kernel::Quintic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_1D * h_inv

    @fastmath if u < 1.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        u_m13 = ( 1.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 + 15.0 * u_m13*u_m13*u_m13*u_m13*u_m13 ) * n
    elseif u < 2.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 ) * n
    elseif u < 1.0
        u_m1  = ( 1.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_1D(kernel::Quintic, u::Real, h_inv::Real)

Evaluate the derivative of the Quintic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_1D(kernel::Quintic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_1D * h_inv^2

    @fastmath if u < 1.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        u_m13 = ( 1.0/3.0 - u )
        return ( -5.0 * u_m1*u_m1*u_m1*u_m1 + 30.0 * u_m23*u_m23*u_m23*u_m23  - 75.0 * u_m13*u_m13*u_m13*u_m13 ) * n
    elseif u < 2.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        return ( -5.0 * u_m1*u_m1*u_m1*u_m1 + 30.0 * u_m23*u_m23*u_m23*u_m23 - 75.0 ) * n
    elseif u < 1.0
        u_m1  = ( 1.0 - u )
        return ( -5.0 * u_m1*u_m1*u_m1*u_m1 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_1D(kernel::Quintic, density::Real, m::Real, h_inv::Real)

Does not do anything for the BSplines. Implemented for stability.
"""
@inline function bias_correction_1D(kernel::Quintic, density::Real, m::Real, h_inv::Real)  
    return density
end

"""
    kernel_value_2D(kernel::Quintic, u::Real, h_inv::Real)

Evaluate quintic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_2D(kernel::Quintic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_2D * h_inv^2

    @fastmath if u < 1.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        u_m13 = ( 1.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 + 15.0 * u_m13*u_m13*u_m13*u_m13*u_m13 ) * n
    elseif u < 2.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 ) * n
    elseif u < 1.0
        u_m1  = ( 1.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_2D(kernel::Quintic, u::Real, h_inv::Real)

Evaluate the derivative of the Quintic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_2D(kernel::Quintic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_2D * h_inv^3

    @fastmath if u < 1.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        u_m13 = ( 1.0/3.0 - u )
        return ( -5.0 * u_m1*u_m1*u_m1*u_m1 + 30.0 * u_m23*u_m23*u_m23*u_m23  - 75.0 * u_m13*u_m13*u_m13*u_m13 ) * n
    elseif u < 2.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        return ( -5.0 * u_m1*u_m1*u_m1*u_m1 + 30.0 * u_m23*u_m23*u_m23*u_m23 - 75.0 ) * n
    elseif u < 1.0
        u_m1  = ( 1.0 - u )
        return ( -5.0 * u_m1*u_m1*u_m1*u_m1 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_2D(kernel::Quintic, density::Real, m::Real, h_inv::Real)

Does not do anything for the BSplines. Implemented for stability.
"""
@inline function bias_correction_2D(kernel::Quintic, density::Real, m::Real, h_inv::Real)  
    return density
end

"""
    kernel_value_3D(kernel::Quintic, u::Real, h_inv::Real)

Evaluate quintic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_3D(kernel::Quintic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_3D * h_inv^3

    @fastmath if u < 1.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        u_m13 = ( 1.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 + 15.0 * u_m13*u_m13*u_m13*u_m13*u_m13 ) * n
    elseif u < 2.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 ) * n
    elseif u < 1.0
        u_m1  = ( 1.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_3D(kernel::Quintic, u::Real, h_inv::Real)

Evaluate the derivative of the Quintic spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_3D(kernel::Quintic, u::Real, h_inv::Real)

    @fastmath n = kernel.norm_3D * h_inv^4

    @fastmath if u < 1.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        u_m13 = ( 1.0/3.0 - u )
        return ( -5.0 * u_m1*u_m1*u_m1*u_m1 + 30.0 * u_m23*u_m23*u_m23*u_m23  - 75.0 * u_m13*u_m13*u_m13*u_m13 ) * n
    elseif u < 2.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        return ( -5.0 * u_m1*u_m1*u_m1*u_m1 + 30.0 * u_m23*u_m23*u_m23*u_m23 - 75.0 ) * n
    elseif u < 1.0
        u_m1  = ( 1.0 - u )
        return ( -5.0 * u_m1*u_m1*u_m1*u_m1 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_3D(kernel::Quintic, density::Real, m::Real, h_inv::Real)

Does not do anything for the BSplines. Implemented for stability.
"""
@inline function bias_correction_3D(kernel::Quintic, density::Real, m::Real, h_inv::Real)  
    return density
end