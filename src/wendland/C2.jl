
"""
    struct WendlandC2{T} <: SPHKernel
        n_neighbours::Int64
        norm_1D::T
        norm_2D::T
        norm_3D::T
    end
"""
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
WendlandC2(T::DataType=Float64, n_neighbours::Integer=100) = WendlandC2{T}(n_neighbours, T(5.0/4.0), T(7.0/π), T(21.0/(2π)))

"""
    kernel_value_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

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
@inline function kernel_deriv_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

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
@inline function bias_correction_1D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T
    n = kernel.norm_1D * h_inv^3
    wc_correction = T(0.0294) * ( kernel.n_neighbours * T(0.01) )^T(-0.977) * m * n
    return density - wc_correction |> T
end

"""
    kernel_value_2D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_2D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

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
@inline function kernel_deriv_2D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

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
@inline function bias_correction_2D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T
    n = kernel.norm_2D * h_inv^3
    wc_correction = T(0.0294) * ( kernel.n_neighbours * T(0.01) )^T(-0.977) * m * n
    return density - wc_correction |> T
end


"""
    kernel_value_3D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
@inline function kernel_value_3D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

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
@inline function kernel_deriv_3D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T

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
@inline function bias_correction_3D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T

    n = kernel.norm_3D * h_inv^3
    wc_correction = T(0.0294) * ( kernel.n_neighbours * T(0.01) )^T(-0.977) * m * n
    return density - wc_correction |> T
end
