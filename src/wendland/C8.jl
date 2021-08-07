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