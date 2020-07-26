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
struct MyKernel <: SPHKernel
    n_neighbours::Int64
    norm_2D::Float64
    norm_3D::Float64
    function MyKernel(n_neighbours::Int64=1000)
        new(n_neighbours, 1.0, 1.0)
    end
end
```

and defining its value and derivative in 2D and 3D, e.g.

```@example 1
@inline function kernel_value_2D(kernel::MyKernel, u::Float64, h_inv::Float64)

    if u < 1.0
        n = kernel.norm_2D * h_inv^2
        return 1.0 * n
    else
        return 0.0
    end

end
```

```@example 1
@inline function kernel_deriv_2D(kernel::MyKernel, u::Float64, h_inv::Float64)
    return 0.0
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