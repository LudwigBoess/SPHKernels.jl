# Adding Kernels

```@meta
CurrentModule = SPHKernels
DocTestSetup = quote
    using SPHKernels
end
```

If you need a different kernel function than the ones I implemented you can add them by defining a new kernel `struct` as a subtype of [SPHKernel](@ref)

```@example 1
using SPHKernels # hide
struct MyKernel{T} <: SPHKernel
    n_neighbours::Int64
    norm_1D::T
    norm_2D::T
    norm_3D::T
end

"""
    MyKernel(T::DataType=Float64, n_neighbours::Integer=295)

Set up a `MyKernel` kernel for a given DataType `T`.
"""
MyKernel(T::DataType=Float64, n_neighbours::Integer=1000) = MyKernel{T}(n_neighbours, T(1), T(1), T(1))
```

and defining its value and derivative in 2D and 3D, e.g.

```@example 1
function kernel_value_2D(kernel::MyKernel{T}, u::Real, h_inv::Real) where T

    if u < 1
        n = kernel.norm_2D * h_inv^2
        return 1.0 * n |> T
    else
        return 0.0 |> T
    end

end
```

```@example 1
@inline function kernel_deriv_2D(kernel::MyKernel{T}, u::Real, h_inv::Real) where T
    return 0.0 |> T
end
```

so you can run the following code

```@example 1
k = MyKernel()
u = 0.5
h_inv = 1.0

v = kernel_value_2D(k, u, h_inv)
println("MyKernel value: $v")

d = kernel_deriv_2D(k, u, h_inv)
println("MyKernel derivative: $d")
```

## Contributing

Please feel free to create a pull request if you feel that your kernel could be useful to others! I only ask you to also add unit tests for your kernels in the `test/runtests.jl` file.