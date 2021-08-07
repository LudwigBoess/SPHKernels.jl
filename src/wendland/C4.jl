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
