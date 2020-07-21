
"""
    struct WendlandC4 <: SPHKernel
        n_neighbours::Integer
        norm_2D::Real
        norm_3D::Real
        function WendlandC4(n_neighbours::Integer=216)
            new(n_neighbours,9.0/π, 495.0/(32.0 * π))
        end
    end
"""
struct WendlandC4 <: SPHKernel
    n_neighbours::Integer
    norm_2D::Real
    norm_3D::Real
    function WendlandC4(n_neighbours::Integer=216)
        new(n_neighbours, 9.0/π, 495.0/(32.0 * π))
    end
end

"""
    kernel_value_2D(kernel::WendlandC4, u::Real, h_inv::Real)

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_2D(kernel::WendlandC4, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^2
        u_m1 = 1.0 - u
        u_m1_5 = u_m1*u_m1*u_m1*u_m1*u_m1
        return ( u_m1_5 * ( 1.0 + 6u + 35.0/3.0 * u^2 ) ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC4, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_2D(kernel::WendlandC4, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^3
        u_m1 = 1.0 - u
        u_m1_5 = u_m1*u_m1*u_m1*u_m1*u_m1
        return ( -288.0/3.0 * u_m1_5 * u^2 - 56.0/3.0 * u * u_m1_5 ) * n
    else
        return 0.
    end

end

"""
    kernel_value_3D(kernel::WendlandC4, u::Real, h_inv::Real)

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_3D(kernel::WendlandC4, u::Real, h_inv::Real)

    @fastmath  if u < 1.0
        n = kernel.norm_3D * h_inv^3
        u_m1 = 1.0 - u
        u_m1_5 = u_m1*u_m1*u_m1*u_m1*u_m1
        return ( u_m1_5 * ( 1.0 + 6u + 35.0/3.0 * u^2 ) ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_3D(kernel::WendlandC4, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_3D(kernel::WendlandC4, u::Real, h_inv::Real)

    if u < 1.0
        n = kernel.norm_3D * h_inv^4
        u_m1 = 1.0 - u
        u_m1_5 = u_m1*u_m1*u_m1*u_m1*u_m1
        return ( -288.0/3.0 * u_m1_5 * u^2 - 56.0/3.0 * u * u_m1_5 ) * n
    else
        return 0.
    end

end


"""
    struct WendlandC6 <: SPHKernel
        n_neighbours::Integer
        norm_2D::Real
        norm_3D::Real
        function WendlandC6(n_neighbours::Integer=295)
            new(n_neighbours, 78.0/(7.0*π), 1365.0/(64.0*π))
        end
    end
"""
struct WendlandC6 <: SPHKernel
    n_neighbours::Integer
    norm_2D::Real
    norm_3D::Real
    function WendlandC6(n_neighbours::Integer=295)
        new(n_neighbours, 78.0/(7.0*π), 1365.0/(64.0*π))
    end
end

"""
    kernel_value_2D(kernel::WendlandC6, u::Real, h_inv::Real)

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_2D(kernel::WendlandC6, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^2
        u_m1 = (1.0 - u)
        u_m1 = u_m1 * u_m1  # (1.0 - u)^2
        u_m1 = u_m1 * u_m1  # (1.0 - u)^4
        u_m1 = u_m1 * u_m1  # (1.0 - u)^8
        u2 = u*u
        return ( u_m1 * ( 1.0 + 8u + 25u2 + 32u2*u )) * n
    else
        return 0.0
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC6, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_2D(kernel::WendlandC6, u::Real, h_inv::Real)


    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^3
        u_m1 = 1.0 - u
        u_m1_7 = u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1
        return ( -22u_m1_7 * u * ( 16u^2 + 7u + 1.0 )) * n
    else
        return 0.0
    end

end

"""
    kernel_value_3D(kernel::WendlandC6, u::Real, h_inv::Real)

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_3D(kernel::WendlandC6, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_3D * h_inv^3
        u_m1 = 1.0 - u
        u_m1 = u_m1 * u_m1  # (1.0 - u)^2
        u_m1 = u_m1 * u_m1  # (1.0 - u)^4
        u_m1 = u_m1 * u_m1  # (1.0 - u)^8
        u2 = u*u
        return ( u_m1 * ( 1.0 + 8u + 25u2 + 32u2*u )) * n
    else
        return 0.0
    end

end

"""
    kernel_deriv_3D(kernel::WendlandC6, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_3D(kernel::WendlandC6, u::Real, h_inv::Real)


    @fastmath if u < 1.0
        n = kernel.norm_3D * h_inv^4
        u_m1 = 1.0 - u
        u_m1_7 = u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1
        return ( -22u_m1_7 * u * ( 16u^2 + 7u + 1.0 )) * n
    else
        return 0.0
    end

end