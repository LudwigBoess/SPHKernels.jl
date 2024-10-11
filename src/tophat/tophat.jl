struct Tophat{T} <: AbstractSPHKernel
    value::T
end

"""
    Tophat(dim::Integer)

Define `Tophat` kernel with dimension `dim` for the native `DataType` of the OS.
"""
Tophat(value::Real=1.0) = Tophat{typeof(value)}(value)

"""
    kernel_value(kernel::Tophat{T}, u::Real) where T

Evaluate Tophat spline at position ``u = \\frac{x}{h}``, without normalisation.
"""
function kernel_value(kernel::Tophat{T}, u::Real) where {T}
    if u < 1
        return kernel.value
    else
        return 0 |> T
    end
end

"""
    kernel_value(kernel::Tophat{T}, u::Real, h_inv::Real) where T

Evaluate Tophat spline at position ``u = \\frac{x}{h}``.
"""
kernel_value(kernel::Tophat{T}, u::Real, h_inv::Real) where {T} = 
    kernel_value(kernel, u)

"""
    kernel_deriv_1D(kernel::Tophat{T}, u::Real) where T

Evaluate the derivative of the Tophat spline at position ``u = \\frac{x}{h}``, without normalisation.
"""
function kernel_deriv(kernel::Tophat{T}, u::Real) where {T}
    return 0 |> T 
end


"""
    kernel_deriv_1D(kernel::Tophat{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the Tophat spline at position ``u = \\frac{x}{h}``.
"""
kernel_deriv(kernel::Tophat{T}, u::Real, h_inv::Real) where {T} = 
    kernel_deriv(kernel, u)

""" 
    bias_correction(kernel::Tophat{T},
                    density::Real, m::Real, h_inv::Real,
                    n_neighbours::Integer) where {T}

Corrects the density estimate for the kernel bias. Not implemented for [`Tophat`](@ref).
"""
function bias_correction(kernel::Tophat{T},
                         density::Real, m::Real, h_inv::Real,
                         n_neighbours::Integer) where {T}
    return density |> T
end