struct WendlandC2{T} <: SPHKernel
    n_neighbours::Int64
    norm_1D::T
    norm_2D::T
    norm_3D::T
end

"""
    WendlandC2(T::DataType=Float64, n_neighbours::Integer=100)

Set up a `WendlandC2` with a given DataType `T`.
"""
WendlandC2(T::DataType=Float64, n_neighbours::Integer=100) = WendlandC2{T}(n_neighbours, 5.0/4.0, 7.0/π, 21.0/(2π))

"""
    kernel_value_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_1D * h_inv
        t1 = 1 - u
        t3 = t1*t1*t1
        return ( t3 * ( 1 + 3u )) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_1D * h_inv^2
        t1 = 1 - u
        return ( -12u * t1  ) * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_1D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_1D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T
    n = kernel.norm_1D * h_inv^3
    wc_correction = T(0.0294) * ( kernel.n_neighbours * T(0.01) )^T(-0.977) * m * n
    return density - wc_correction |> T
end

"""
    kernel_value_2D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_2D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^2
        t1 = 1 - u
        t4 = t1 * t1 * t1 * t1
        return ( t4 * ( 1 + 4u ) ) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_2D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^3
        t1 = 1 - u
        t3 = t1 * t1 * t1
        return ( -20u * t3 ) * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_2D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_2D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T
    n = kernel.norm_2D * h_inv^3
    wc_correction = T(0.0294) * ( kernel.n_neighbours * T(0.01) )^T(-0.977) * m * n
    return density - wc_correction |> T
end


"""
    kernel_value_3D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_3D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

    if u < 1
        n  = kernel.norm_3D * h_inv^3
        t1 = 1 - u
        t4 = t1 * t1 * t1 * t1
        return ( t4 * ( 1 + 4u ) ) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_3D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_3D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_3D * h_inv^2 * h_inv^2
        t1 = 1 - u
        t3 = t1 * t1 * t1
        return ( -20u * t3 ) * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_3D(kernel::WendlandC2, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_3D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T

    n = kernel.norm_3D * h_inv^3
    wc_correction = T(0.0294) * ( kernel.n_neighbours * T(0.01) )^T(-0.977) * m * n
    return density - wc_correction |> T
end



struct WendlandC4{T} <: SPHKernel
    n_neighbours::Int64
    norm_1D::T
    norm_2D::T
    norm_3D::T
end

"""
    WendlandC4(T::DataType=Float64, n_neighbours::Integer=216)

Set up a `WendlandC4` kernel for a given DataType `T`.
"""
WendlandC4(T::DataType=Float64, n_neighbours::Integer=216) = WendlandC4{T}(n_neighbours, 3.0/2.0, 9.0/π, 495.0/(32.0 * π))

"""
    kernel_value_1D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_1D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_1D * h_inv
        t1 = 1 - u
        t5 = t1*t1*t1*t1*t1
        return ( t5 * (1 + 5u + 8u^2 )) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_1D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_1D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_1D * h_inv^2
        t1 = 1 - u
        t4 = t1*t1*t1*t1
        return ( -14u * t4 - 56u^2 * t4 ) * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_1D(kernel::WendlandC4{T}, density::Real, m::Real, h_inv::Real) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_1D(kernel::WendlandC4{T}, density::Real, m::Real, h_inv::Real) where T
    n = kernel.norm_1D * h_inv^3
    wc_correction = T(0.01342) * ( kernel.n_neighbours * T(0.01) )^T(-1.579) * m * n
    
    return density - wc_correction |> T
end

"""
    (kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_2D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^2
        t1 = 1 - u
        t5 = t1*t1*t1*t1*t1
        return ( t1 * t5 * ( 1 + 6u + 35/3 * u^2 ) ) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_2D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^3
        t1 = 1 - u
        t5 = t1*t1*t1*t1*t1
        return ( -288/3 * t5 * u^2 - 56/3 * u * t5 ) * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_2D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_2D(kernel::WendlandC4{T}, density::Real, m::Real, h_inv::Real) where T

    n = kernel.norm_2D * h_inv^3
    wc_correction = T(0.01342) * ( kernel.n_neighbours * T(0.01) )^T(-1.579) * m * n
    return density - wc_correction
end


"""
    kernel_value_3D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_3D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_3D * h_inv^3
        t1 = 1 - u
        t5 = t1*t1*t1*t1*t1
        return ( t1 * t5 * ( 1 + 6u + 35/3 * u^2 ) ) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_3D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC4 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_3D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_3D * h_inv^2 * h_inv^2
        t1 = 1 - u
        t5 = t1*t1*t1*t1*t1
        return ( -288/3 * t5 * u^2 - 56/3 * u * t5 ) * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_3D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_3D(kernel::WendlandC4{T}, density::Real, m::Real, h_inv::Real) where T

    n = kernel.norm_3D * h_inv^3
    wc_correction = T(0.01342) * ( kernel.n_neighbours * T(0.01) )^T(-1.579) * m * n
    
    return density - wc_correction |> T
end



struct WendlandC6{T} <: SPHKernel
    n_neighbours::Int64
    norm_1D::T
    norm_2D::T
    norm_3D::T
end

"""
    WendlandC6(T::DataType=Float64, n_neighbours::Integer=295)

Set up a `WendlandC6` kernel for a given DataType `T`.
"""
WendlandC6(T::DataType=Float64, n_neighbours::Integer=295) = WendlandC6{T}(n_neighbours, T(55.0/32.0), T(78.0/(7.0*π)), T(1365.0/(64.0*π)))


"""
    kernel_value_1D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_1D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_1D * h_inv
        t1 = 1 - u
        t6 = t1*t1*t1*t1*t1*t1
        u2 = u*u
        return ( t6 * t1 * (1 + 7u + 19u2 + 21u2 * u)) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_1D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_1D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

     if u < 1
        n = kernel.norm_1D * h_inv^2
        t1 = 1 - u
        t6 = t1*t1*t1*t1*t1*t1
        u2 = u*u
        return ( -6t6 * u * (35u2 + 18u + 3 )) * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_1D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_1D(kernel::WendlandC6{T}, density::Real, m::Real, h_inv::Real) where T

    n = kernel.norm_1D * h_inv^3
    wc_correction = T(0.0116) * ( kernel.n_neighbours * T(0.01) )^T(-2.236) * m * n

    if wc_correction < T(0.2)*density
        density -= wc_correction
    end
    
    return density |> T
end


"""
    (kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_2D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^2
        u_m1 = 1 - u
        u_m1 = u_m1 * u_m1  # (1.0 - u)^2
        u_m1 = u_m1 * u_m1  # (1.0 - u)^4
        u_m1 = u_m1 * u_m1  # (1.0 - u)^8
        u2 = u*u
        return ( u_m1 * ( 1 + 8u + 25u2 + 32u2*u )) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC6, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_2D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T


    if u < 1
        n = kernel.norm_2D * h_inv^3
        u_m1 = 1 - u
        u_m1_7 = u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1
        return ( -22u_m1_7 * u * ( 16u^2 + 7u + 1 )) * n |> T
    else
        return 0.0 |> T 
    end

end

""" 
    bias_correction_2D(kernel::WendlandC6{T}, density::Real, m::Real, h_inv::Real) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_2D(kernel::WendlandC6{T}, density::Real, m::Real, h_inv::Real) where T

    n = kernel.norm_2D * h_inv^3
    wc_correction = T(0.0116) * ( kernel.n_neighbours * T(0.01) )^T(-2.236) * m * n

    if wc_correction < T(0.2)*density
        density -= wc_correction
    end
    
    return density |> T
end


"""
    kernel_value_3D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_3D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_3D * h_inv^3
        u_m1 = 1 - u
        u_m1 = u_m1 * u_m1  # (1.0 - u)^2
        u_m1 = u_m1 * u_m1  # (1.0 - u)^4
        u_m1 = u_m1 * u_m1  # (1.0 - u)^8
        u2 = u*u
        return ( u_m1 * ( 1 + 8u + 25u2 + 32u2*u )) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_3D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_3D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_3D * h_inv^2 * h_inv^2
        u_m1 = 1 - u
        u_m1_7 = u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1 * u_m1
        return ( -22u_m1_7 * u * ( 16u^2 + 7u + 1 )) * n |> T
    else
        return 0.0 |> T
    end

end


""" 
    bias_correction_3D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_3D(kernel::WendlandC6{T}, density::Real, m::Real, h_inv::Real) where T

    n = kernel.norm_3D * h_inv^3
    wc_correction::T = T(0.0116) * ( kernel.n_neighbours * T(0.01) )^T(-2.236) * m * n

    if wc_correction < T(0.2)*density
        density -= wc_correction
    end
    
    return density |> T
end




struct WendlandC8{T} <: SPHKernel
    n_neighbours::Int64
    norm_2D::T
    norm_3D::T 
end

"""
    WendlandC8(T::DataType=Float64, n_neighbours::Integer=395)

Set up a `WendlandC8` kernel for a given DataType `T`.
"""
WendlandC8(T::DataType=Float64, n_neighbours::Integer=395) = WendlandC8{T}(n_neighbours, 8.0/(3π), 357.0/(64π))
    

"""
    kernel_value_2D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_2D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^2
        t1 = 1 - u
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return ( t1 * t9 * (5 + 50u + 210u2 + 450u2 * u + 429u2 * u2)) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_2D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC6 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_2D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^3
        t1 = 1.0 - u
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return -26t9 * u * (231u2 * u + 159u2 + 45u + 5 ) * n |> T
    else
        return 0.0 |> T
    end

end

""" 
    bias_correction_2D(kernel::WendlandC8, density::Real, m::Real, h_inv::Real)

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_2D(kernel::WendlandC8{T}, density::Real, m::Real, h_inv::Real) where T
    return density |> T
end


"""
    kernel_value_3D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC8 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value_3D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^3
        t1 = (1 - u)
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return ( t1 * t9 * (5 + 50u + 210u2 + 450u2 * u + 429u2 * u2)) * n |> T
    else
        return 0.0 |> T
    end

end

"""
    kernel_deriv_3D(kernel::WendlandC8, u::Real, h_inv::Real)

Evaluate the derivative of the WendlandC8 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv_3D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^2 * h_inv^2
        t1 = 1 - u
        t9 = t1*t1*t1*t1*t1*t1*t1*t1*t1  # (1.0 - u)^9
        u2 = u*u
        return -26t9 * u * (231u2 * u + 159u2 + 45u + 5 ) * n |> T
    else
        return 0.0 |> T
    end

end


""" 
    bias_correction_3D(kernel::WendlandC8{T}, density::Real, m::Real, h_inv::Real) where T

Corrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
"""
function bias_correction_3D(kernel::WendlandC8{T}, density::Real, m::Real, h_inv::Real) where T
    return density |> T
end