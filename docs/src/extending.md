# Adding Kernels

```@meta
CurrentModule = SPHKernels
DocTestSetup = quote
    using SPHKernels
end
```

If you need a different kernel function than the ones I implemented you can add them by defining a new kernel `struct` as a subtype of [AbstractSPHKernel](@ref)

```@example 1
using SPHKernels # hide
struct MyKernel{T} <: AbstractSPHKernel
    dim::Int64
    norm::T
end

"""
    MyKernel(T::DataType=Float64, dim::Integer=3)

Set up a `MyKernel` kernel for a given DataType `T`.
"""
MyKernel(T::DataType=Float64, dim::Integer=3) = MyKernel{T}(dim, T(1))


"""
    MyKernel(dim::Integer)

Define `MyKernel` kernel with dimension `dim` for the native `DataType` of the OS.
"""
MyKernel(dim::Integer) = MyKernel{T}(typeof(1.0), dim)
```

and defining its value and derivative, e.g.

```@example 1
function kernel_value(kernel::MyKernel{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm * h_inv^kernel.dim
        return 1.0 * n |> T
    else
        return 0.0 |> T
    end

end
```

```@example 1
@inline function kernel_deriv(kernel::MyKernel{T}, u::Real, h_inv::Real) where T
    return 0.0 |> T
end
```

so you can run the following code

```@example 1
k = MyKernel()
u = 0.5
h_inv = 1.0

v = kernel_value(k, u, h_inv)
println("MyKernel value: $v")

d = kernel_derivk, u, h_inv)
println("MyKernel derivative: $d")
```

## Contributing

Please feel free to create a pull request if you feel that your kernel could be useful to others! I only ask you to also add unit tests for your kernels in the `test/runtests.jl` file.