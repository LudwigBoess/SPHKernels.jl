
"""
    struct WendlandC2 <: SPHKernel
        n_neighbours::Int64
        norm_1D::Float64
        norm_2D::Float64
        norm_3D::Float64
    end
"""
struct WendlandC2 <: SPHKernel
    n_neighbours::Int64
    norm_1D::Float64
    norm_2D::Float64
    norm_3D::Float64
    function WendlandC2(n_neighbours::Integer=100)
        new(n_neighbours, 5.0/4.0, 7.0/π, 21.0/(2π))
    end
end

"""
    kernel_value_1D(kernel::WendlandC2, u::Real, h_inv::Real)

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_1D(kernel::WendlandC2, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_1D * h_inv
        t1 = 1.0 - u
        t3 = t1*t1*t1
        return ( t3 * (1.0 + 3u )) * n
    else
        return 0.0
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC2, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_1D(kernel::WendlandC2, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_1D * h_inv^2
        t1 = 1.0 - u
        return ( -12.0 * u * t  ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_1D(kernel::WendlandC2, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_1D(kernel::WendlandC2, density::Real, m::Real, h_inv::Real)
    return density
end

"""
    kernel_value_2D(kernel::WendlandC2, u::Real, h_inv::Real)

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_2D(kernel::WendlandC2, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^2
        t1 = 1.0 - u
        t4 = t * t * t * t
        return ( t4 * ( 1.0 + 4u ) ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC2, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_2D(kernel::WendlandC2, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^3
        t1 = 1.0 - u
        t3 = t1 * t1 * t1
        return ( -20.0 * u * t3 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_2D(kernel::WendlandC2, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_2D(kernel::WendlandC2, density::Real, m::Real, h_inv::Real)
    return density 
end


"""
    kernel_value_3D(kernel::WendlandC2, u::Real, h_inv::Real)

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_3D(kernel::WendlandC2, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^3
        t1 = 1.0 - u
        t4 = t * t * t * t
        return ( t4 * ( 1.0 + 4u ) ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_3D(kernel::WendlandC2, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_3D(kernel::WendlandC2, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^4
        t1 = 1.0 - u
        t3 = t1 * t1 * t1
        return ( -20.0 * u * t3 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_3D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_3D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)

    @fastmath n = kernel.norm_3D * h_inv^3
    @fastmath wc_correction = 0.0294 * ( kernel.n_neighbours * 0.01 )^(-1.977) * m * n
    return density - wc_correction
end


"""
    struct WendlandC4 <: SPHKernel
        n_neighbours::Int64
        norm_1D::Float64
        norm_2D::Float64
        norm_3D::Float64
    end
"""
struct WendlandC4 <: SPHKernel
    n_neighbours::Int64
    norm_1D::Float64
    norm_2D::Float64
    norm_3D::Float64
    function WendlandC4(n_neighbours::Integer=216)
        new(n_neighbours, 3.0/2.0, 9.0/π, 495.0/(32.0 * π))
    end
end

"""
    kernel_value_1D(kernel::WendlandC4, u::Real, h_inv::Real)

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_1D(kernel::WendlandC4, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_1D * h_inv
        t1 = 1.0 - u
        t5 = t1*t1*t1*t1*t1
        return ( t5 * (1.0 + 5u + 8u^2 )) * n
    else
        return 0.0
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC4, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_1D(kernel::WendlandC4, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_1D * h_inv^2
        t1 = 1.0 - u
        t4 = t1*t1*t1*t1
        return ( -14.0 * u * t4 - 56.0 * u^2 * t4 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_1D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_1D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)
    return density
end

"""
    kernel_value_2D(kernel::WendlandC4, u::Real, h_inv::Real)

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_2D(kernel::WendlandC4, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^2
        t1 = 1.0 - u
        t5 = t1*t1*t1*t1*t1
        return ( t1 * t5 * ( 1.0 + 6u + 35.0/3.0 * u^2 ) ) * n
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
        t1 = 1.0 - u
        t5 = t1*t1*t1*t1*t1
        return ( -288.0/3.0 * t5 * u^2 - 56.0/3.0 * u * t5 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_2D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_2D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)

    @fastmath n = kernel.norm_2D * h_inv^3
    @fastmath wc_correction = 0.01342 * ( kernel.n_neighbours * 0.01 )^(-1.579) * m * n
    return density - wc_correction
end


"""
    kernel_value_3D(kernel::WendlandC4, u::Real, h_inv::Real)

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_3D(kernel::WendlandC4, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_3D * h_inv^3
        t1 = 1.0 - u
        t5 = t1*t1*t1*t1*t1
        return ( t1 * t5 * ( 1.0 + 6u + 35.0/3.0 * u^2 ) ) * n
    else
        return 0.
    end

end

"""
    kernel_deriv_3D(kernel::WendlandC4, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_3D(kernel::WendlandC4, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_3D * h_inv^4
        t1 = 1.0 - u
        t5 = t1*t1*t1*t1*t1
        return ( -288.0/3.0 * t5 * u^2 - 56.0/3.0 * u * t5 ) * n
    else
        return 0.
    end

end

""" 
    bias_correction_3D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_3D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)

    @fastmath n = kernel.norm_3D * h_inv^3
    @fastmath wc_correction = 0.01342 * ( kernel.n_neighbours * 0.01 )^(-1.579) * m * n
    
    return density - wc_correction
end


"""
    struct WendlandC6 <: SPHKernel
        n_neighbours::Int64
        norm_1D::Float64
        norm_2D::Float64
        norm_3D::Float64
    end
"""
struct WendlandC6 <: SPHKernel
    n_neighbours::Int64
    norm_1D::Float64
    norm_2D::Float64
    norm_3D::Float64
    function WendlandC6(n_neighbours::Integer=295)
        new(n_neighbours, 55.0/32.0, 78.0/(7.0*π), 1365.0/(64.0*π))
    end
end


"""
    kernel_value_1D(kernel::WendlandC6, u::Real, h_inv::Real)

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_1D(kernel::WendlandC6, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_1D * h_inv
        t1 = 1.0 - u
        t6 = t1*t1*t1*t1*t1*t1
        u2 = u*u
        return ( t6 * t1 * (1.0 + 7u + 19u2 + 21u2 * u)) * n
    else
        return 0.0
    end

end

"""
    kernel_deriv_1D(kernel::WendlandC6, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_1D(kernel::WendlandC6, u::Real, h_inv::Real)


    @fastmath if u < 1.0
        n = kernel.norm_1D * h_inv^2
        t1 = 1.0 - u
        t6 = t1*t1*t1*t1*t1*t1
        u2 = u*u
        return ( -6t6 * u * (35u2 + 18u + 3.0)) * n
    else
        return 0.0
    end

end

""" 
    bias_correction_1D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_1D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)

    @fastmath n = kernel.norm_1D * h_inv^3
    @fastmath wc_correction = 0.0116 * ( kernel.n_neighbours * 0.01 )^(-2.236) * m * n

    if wc_correction < 0.2*density
        density -= wc_correction
    end
    
    return density
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
    bias_correction_2D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_2D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)

    @fastmath n = kernel.norm_2D * h_inv^3
    @fastmath wc_correction = 0.0116 * ( kernel.n_neighbours * 0.01 )^(-2.236) * m * n

    if wc_correction < 0.2*density
        density -= wc_correction
    end
    
    return density
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


""" 
    bias_correction_3D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_3D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)

    @fastmath n = kernel.norm_3D * h_inv^3
    @fastmath wc_correction = 0.0116 * ( kernel.n_neighbours * 0.01 )^(-2.236) * m * n

    if wc_correction < 0.2*density
        density -= wc_correction
    end
    
    return density
end




"""
    struct WendlandC8 <: SPHKernel
        n_neighbours::Int64
        norm_2D::Float64
        norm_3D::Float64
        function WendlandC6(n_neighbours::Integer=395)
            new(n_neighbours, 8.0/(3π), 357.0/(64π))
        end
    end
"""
struct WendlandC8 <: SPHKernel
    n_neighbours::Int64
    norm_2D::Float64
    norm_3D::Float64
    function WendlandC8(n_neighbours::Integer=395)
        new(n_neighbours, 8.0/(3π), 357.0/(64π))
    end
end

"""
    kernel_value_2D(kernel::WendlandC8, u::Real, h_inv::Real)

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_2D(kernel::WendlandC8, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^2
        t1 = (1.0 - u)
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return ( t1 * t9 * (5.0 + 50u + 210u2 + 450u2 * u + 429u2 * u2)) * n
    else
        return 0.0
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC8, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_2D(kernel::WendlandC8, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^2
        t1 = (1.0 - u)
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return -26t9 * u * (231u2 * u + 159u2 + 45u + 5.0) * n
    else
        return 0.0
    end

end

""" 
    bias_correction_2D(kernel::WendlandC8, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_2D(kernel::WendlandC8, density::Real, m::Real, h_inv::Real)

    return density
end


"""
    kernel_value_3D(kernel::WendlandC8, u::Real, h_inv::Real)

Evaluate WendlandC8 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_3D(kernel::WendlandC8, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^2
        t1 = (1.0 - u)
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return ( t1 * t9 * (5.0 + 50u + 210u2 + 450u2 * u + 429u2 * u2)) * n
    else
        return 0.0
    end

end

"""
    kernel_deriv_3D(kernel::WendlandC8, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC8 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_deriv_3D(kernel::WendlandC8, u::Real, h_inv::Real)

    @fastmath if u < 1.0
        n = kernel.norm_2D * h_inv^2
        t1 = (1.0 - u)
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return -26t9 * u * (231u2 * u + 159u2 + 45u + 5.0) * n
    else
        return 0.0
    end

end


""" 
    bias_correction_3D(kernel::WendlandC8, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
@inline function bias_correction_3D(kernel::WendlandC8, density::Real, m::Real, h_inv::Real)
    return density
end
