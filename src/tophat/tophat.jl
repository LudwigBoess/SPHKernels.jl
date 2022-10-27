struct Tophat{T} <: AbstractSPHKernel
    value::T
end


"""
    WendlandC2(dim::Integer)

Define `WendlandC2` kernel with dimension `dim` for the native `DataType` of the OS.
"""
Tophat(value::Real=1.0) = Tophat{typeof(value)}(value)

"""
    kernel_value(kernel::Tophat{T}, u::Real, h_inv::Real) where T

Evaluate WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::Tophat{T}, u::Real, h_inv::Real) where {T}

    if u < 1
        return kernel.value
    else
        return 0 |> T
    end
end

"""
    kernel_deriv_1D(kernel::Tophat{T}, u::Real, h_inv::Real) where T

Evaluate the derivative of the WendlandC2 spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::Tophat{T}, u::Real, h_inv::Real) where {T}
    return 0 |> T 
end


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