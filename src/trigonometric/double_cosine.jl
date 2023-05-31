"""
    Double Cosine kernel as in Yang+13, https://doi.org/10.1016/j.apm.2013.12.001
"""

"""
    struct DoubleCosine{T} <: AbstractSPHKernel
        dim::Int64
        norm::T
    end
"""
struct DoubleCosine{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end


"""
    DoubleCosine(T::DataType=Float64, dim::Integer=3)

Set up a `DoubleCosine` kernel for a given DataType `T`.
"""
function DoubleCosine(T::DataType=Float64, dim::Integer=3)

    if dim == 1
        norm = T(1 / 6)
    elseif dim == 2
        norm = T(π / (3π^2 - 16))
    elseif dim == 3
        norm = T(π^2 / (4π^2 - 30))
    else
        error("DoubleCosine not defined for $dim dimensions!")
    end

    return DoubleCosine{T}(dim, norm)
end

"""
    DoubleCosine(dim::Integer)

Define `DoubleCosine` kernel with dimension `dim` for the native `DataType` of the OS.
"""
DoubleCosine(dim::Integer) = DoubleCosine(typeof(1.0), dim)

"""
    kernel_value(kernel::DoubleCosine, u::Real, h_inv::Real)
Evaluate DoubleCosine spline at position ``u = \\frac{x}{h}``.
"""
function kernel_value(kernel::DoubleCosine{T}, u::Real, h_inv::Real) where {T}

    if u < 1
        n = kernel.norm * h_inv^kernel.dim
        return (4 * cos(π * u) + cos(2π * u) + 3) * n
    else
        return 0 |> T
    end

end


"""
    kernel_deriv(kernel::DoubleCosine, u::Real, h_inv::Real)

Evaluate the derivative of the DoubleCosine spline at position ``u = \\frac{x}{h}``.
"""
function kernel_deriv(kernel::DoubleCosine{T}, u::Real, h_inv::Real) where {T}

    if u < 1
        n = kernel.norm * h_inv^kernel.dim * h_inv
        return (-4π * sin(π * u) - 2π * sin(2π * u)) * n
    else
        return 0 |> T
    end

end


""" 
    bias_correction( kernel::DoubleCosine{T}, 
                     density::Real, m::Real, h_inv::Real, 
                     n_neighbours::Integer )  where T

Does not do anything for the DoubleCosine. Implemented for stability.
"""
function bias_correction(kernel::DoubleCosine{T},
                        density::Real, m::Real, h_inv::Real,
                        n_neighbours::Integer) where {T}
    return density |> T
end