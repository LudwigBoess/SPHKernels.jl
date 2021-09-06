var documenterSearchIndex = {"docs":
[{"location":"kernels/#Kernels","page":"Kernels","title":"Kernels","text":"","category":"section"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"CurrentModule = SPHKernels\nDocTestSetup = quote\n    using SPHKernels\nend","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"This package supplies a number of kernels frequently used in Smoothed-Particle Hydrodynamics (SPH), as well as functions to evaluate their values and derivatives in 2D and 3D.","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"These kernels include the B-splines (Cubic and Quintic) suggested in Monaghan & Lattanzio (1985) and the Wendland functions (WendlandC2, WendlandC4 and WendlandC6) as suggested in Dehnen & Aly (2012).","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"In this implementation we follow the convention of Dehnen&Aly in using the 'compact kernel support' as a means to define the maximum extent of the kernel. They denote this H in their paper, for convenience (aka for not having to type caps) we use h in the code.","category":"page"},{"location":"kernels/#Evaluating-Kernels","page":"Kernels","title":"Evaluating Kernels","text":"","category":"section"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"To evaluate a 3D kernel you need to use the function","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"kernel_value_3D(k::SPHKernel, u::Real, h_inv::Real)","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"where SPHKernel is the supertype for an implemented SPH kernel, u = fracxh is the distance to the kernel origin in measures of the compact kernel support and h_inv is the inverse of the compact kernel support.","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"The same goes for a 1D kernel","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"kernel_value_1D(k::SPHKernel, u::Real, h_inv::Real)","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"and a 2D kernel","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"kernel_value_2D(k::SPHKernel, u::Real, h_inv::Real)","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"If you want your code to look a little more fancy you can also use the alternative functions 𝒲₁, where the respective subscript refers to the dimension:","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"𝒲₁( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value_1D(kernel, u, h_inv)\n𝒲₂( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value_2D(kernel, u, h_inv)\n𝒲₃( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value_3D(kernel, u, h_inv)","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"As an example:","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"using SPHKernels # hide\nk     = WendlandC6()\n# distance between the particle and the origin of the kernel\nr     = 0.5\nh     = 1.0\nh_inv = 1.0/h\nu     = r * h_inv\n\n# kernel value at position r\nval = 𝒲₃(k, u, h_inv)\n\nprintln(\"val = $val\")","category":"page"},{"location":"kernels/#Evaluating-Derivatives","page":"Kernels","title":"Evaluating Derivatives","text":"","category":"section"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"Similar to Evaluating Kernels you can evluate a kernel derivative with","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"kernel_deriv_3D(k::SPHKernel, u::Real, h_inv::Real)","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"or in the fancy way:","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"∇𝒲₁(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv_1D(kernel, u, h_inv)\n∇𝒲₂(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv_2D(kernel, u, h_inv)\n∇𝒲₃(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv_3D(kernel, u, h_inv)\n","category":"page"},{"location":"kernels/#Bias-Correction","page":"Kernels","title":"Bias Correction","text":"","category":"section"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"You can correct for the kernel bias of the Wendland kernels as described in Dehnen & Aly (2012), Eq. 18 + 19 with the functions:","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"bias_correction_1D(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)\nbias_correction_2D(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)\nbias_correction_3D(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"or again in the fancy way","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"δρ₁(kernel::SPHKernel, density::Real, m::Real, h_inv::Real) = bias_correction_1D(kernel, density, m, h_inv)\nδρ₂(kernel::SPHKernel, density::Real, m::Real, h_inv::Real) = bias_correction_2D(kernel, density, m, h_inv)\nδρ₃(kernel::SPHKernel, density::Real, m::Real, h_inv::Real) = bias_correction_3D(kernel, density, m, h_inv)\n","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"This will return a new value for the density:","category":"page"},{"location":"kernels/","page":"Kernels","title":"Kernels","text":"using SPHKernels # hide\ndensity = 1.0\nkernel  = WendlandC6()\n\n# correct density\ndensity = bias_correction_3D(kernel, density, 1.0, 0.5)\n\nprintln(\"density = $density\")","category":"page"},{"location":"api/#API-Reference","page":"API reference","title":"API Reference","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"CurrentModule = SPHKernels\nDocTestSetup = quote\n    using SPHKernels\nend","category":"page"},{"location":"api/","page":"API reference","title":"API reference","text":"","category":"page"},{"location":"api/#Exported-Types","page":"API reference","title":"Exported Types","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [SPHKernels]\nPrivate = false\nOrder = [:type]","category":"page"},{"location":"api/#SPHKernels.Cubic","page":"API reference","title":"SPHKernels.Cubic","text":"Cubic(T::DataType=Float64, n_neighbours::Integer=64)\n\nSet up a Cubic kernel for a given DataType T.\n\n\n\n\n\n","category":"type"},{"location":"api/#SPHKernels.Quintic","page":"API reference","title":"SPHKernels.Quintic","text":"Quintic(T::DataType=Float64, n_neighbours::Integer=64)\n\nSet up a Quintic kernel for a given DataType T.\n\n\n\n\n\n","category":"type"},{"location":"api/#SPHKernels.SPHKernel","page":"API reference","title":"SPHKernels.SPHKernel","text":"SPHKernel\n\nSupertype for all SPH kernels.\n\n\n\n\n\n","category":"type"},{"location":"api/#SPHKernels.WendlandC2","page":"API reference","title":"SPHKernels.WendlandC2","text":"struct WendlandC2{T} <: SPHKernel\n    n_neighbours::Int64\n    norm_1D::T\n    norm_2D::T\n    norm_3D::T\nend\n\n\n\n\n\n","category":"type"},{"location":"api/#SPHKernels.WendlandC2-2","page":"API reference","title":"SPHKernels.WendlandC2","text":"WendlandC2(T::DataType=Float64, n_neighbours::Integer=100)\n\nSet up a WendlandC2 with a given DataType T.\n\n\n\n\n\n","category":"type"},{"location":"api/#SPHKernels.WendlandC4","page":"API reference","title":"SPHKernels.WendlandC4","text":"WendlandC4(T::DataType=Float64, n_neighbours::Integer=216)\n\nSet up a WendlandC4 kernel for a given DataType T.\n\n\n\n\n\n","category":"type"},{"location":"api/#SPHKernels.WendlandC6","page":"API reference","title":"SPHKernels.WendlandC6","text":"WendlandC6(T::DataType=Float64, n_neighbours::Integer=295)\n\nSet up a WendlandC6 kernel for a given DataType T.\n\n\n\n\n\n","category":"type"},{"location":"api/#SPHKernels.WendlandC8","page":"API reference","title":"SPHKernels.WendlandC8","text":"WendlandC8(T::DataType=Float64, n_neighbours::Integer=395)\n\nSet up a WendlandC8 kernel for a given DataType T.\n\n\n\n\n\n","category":"type"},{"location":"api/#Exported-Functions","page":"API reference","title":"Exported Functions","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [SPHKernels]\nPrivate = false\nOrder = [:function]","category":"page"},{"location":"api/#SPHKernels.bias_correction_1D-Union{Tuple{T}, Tuple{Cubic{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_1D","text":"bias_correction_1D(kernel::Cubic{T}, density::Real, m::Real, h_inv::Real) where T\n\nDoes not do anything for the BSplines. Implemented for stability.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_1D-Union{Tuple{T}, Tuple{Quintic{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_1D","text":"bias_correction_1D(kernel::Quintic{T}, density::Real, m::Real, h_inv::Real) where T\n\nDoes not do anything for the BSplines. Implemented for stability.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_1D-Union{Tuple{T}, Tuple{WendlandC2{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_1D","text":"bias_correction_1D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_1D-Union{Tuple{T}, Tuple{WendlandC4{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_1D","text":"bias_correction_1D(kernel::WendlandC4{T}, density::Real, m::Real, h_inv::Real) where T\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_1D-Union{Tuple{T}, Tuple{WendlandC6{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_1D","text":"bias_correction_1D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_2D-Union{Tuple{T}, Tuple{Cubic{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_2D","text":"bias_correction_2D(kernel::Cubic, density::Real, m::Real, h_inv::Real)\n\nDoes not do anything for the BSplines. Implemented for stability.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_2D-Union{Tuple{T}, Tuple{Quintic{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_2D","text":"bias_correction_2D(kernel::Quintic{T}, density::Real, m::Real, h_inv::Real) where T\n\nDoes not do anything for the BSplines. Implemented for stability.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_2D-Union{Tuple{T}, Tuple{WendlandC2{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_2D","text":"bias_correction_2D(kernel::WendlandC2{T}, density::Real, m::Real, h_inv::Real) where T\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_2D-Union{Tuple{T}, Tuple{WendlandC4{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_2D","text":"bias_correction_2D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_2D-Union{Tuple{T}, Tuple{WendlandC6{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_2D","text":"bias_correction_2D(kernel::WendlandC6{T}, density::Real, m::Real, h_inv::Real) where T\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_2D-Union{Tuple{T}, Tuple{WendlandC8{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_2D","text":"bias_correction_2D(kernel::WendlandC8, density::Real, m::Real, h_inv::Real)\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_3D-Union{Tuple{T}, Tuple{Cubic{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_3D","text":"bias_correction_3D(kernel::Cubic{T}, density::Real, m::Real, h_inv::Real) where T\n\nDoes not do anything for the BSplines. Implemented for stability.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_3D-Union{Tuple{T}, Tuple{Quintic{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_3D","text":"bias_correction_3D(kernel::Quintic{T}, density::Real, m::Real, h_inv::Real) where T\n\nDoes not do anything for the BSplines. Implemented for stability.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_3D-Union{Tuple{T}, Tuple{WendlandC2{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_3D","text":"bias_correction_3D(kernel::WendlandC2, density::Real, m::Real, h_inv::Real)\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_3D-Union{Tuple{T}, Tuple{WendlandC4{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_3D","text":"bias_correction_3D(kernel::WendlandC4, density::Real, m::Real, h_inv::Real)\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_3D-Union{Tuple{T}, Tuple{WendlandC6{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_3D","text":"bias_correction_3D(kernel::WendlandC6, density::Real, m::Real, h_inv::Real)\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.bias_correction_3D-Union{Tuple{T}, Tuple{WendlandC8{T}, Real, Real, Real}} where T","page":"API reference","title":"SPHKernels.bias_correction_3D","text":"bias_correction_3D(kernel::WendlandC8{T}, density::Real, m::Real, h_inv::Real) where T\n\nCorrects the density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_1D-Union{Tuple{T}, Tuple{Cubic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_1D","text":"kernel_deriv_1D(kernel::Cubic{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the Cubic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_1D-Union{Tuple{T}, Tuple{Quintic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_1D","text":"kernel_deriv_1D(kernel::Quintic{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the Quintic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_1D-Union{Tuple{T}, Tuple{WendlandC2{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_1D","text":"kernel_deriv_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the WendlandC2 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_1D-Union{Tuple{T}, Tuple{WendlandC4{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_1D","text":"kernel_deriv_1D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the WendlandC4 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_1D-Union{Tuple{T}, Tuple{WendlandC6{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_1D","text":"kernel_deriv_1D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the WendlandC6 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_2D-Union{Tuple{T}, Tuple{Cubic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_2D","text":"kernel_deriv_2D(kernel::Cubic{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the Cubic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_2D-Union{Tuple{T}, Tuple{Quintic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_2D","text":"kernel_deriv_2D(kernel::Quintic{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the Quintic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_2D-Union{Tuple{T}, Tuple{WendlandC2{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_2D","text":"kernel_deriv_2D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the WendlandC2 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_2D-Union{Tuple{T}, Tuple{WendlandC4{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_2D","text":"kernel_deriv_2D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the WendlandC4 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_2D-Union{Tuple{T}, Tuple{WendlandC6{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_2D","text":"kernel_deriv_2D(kernel::WendlandC6, u::Real, h_inv::Real)\n\nEvaluate the derivative of the WendlandC6 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_2D-Union{Tuple{T}, Tuple{WendlandC8{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_2D","text":"kernel_deriv_2D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the WendlandC6 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_3D-Union{Tuple{T}, Tuple{Cubic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_3D","text":"kernel_deriv_3D(kernel::Cubic{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the Cubic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_3D-Union{Tuple{T}, Tuple{Quintic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_3D","text":"kernel_deriv_3D(kernel::Quintic, u::Real, h_inv::Real)\n\nEvaluate the derivative of the Quintic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_3D-Union{Tuple{T}, Tuple{WendlandC2{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_3D","text":"kernel_deriv_3D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the WendlandC2 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_3D-Union{Tuple{T}, Tuple{WendlandC4{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_3D","text":"kernel_deriv_3D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the WendlandC4 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_3D-Union{Tuple{T}, Tuple{WendlandC6{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_3D","text":"kernel_deriv_3D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T\n\nEvaluate the derivative of the WendlandC6 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_deriv_3D-Union{Tuple{T}, Tuple{WendlandC8{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_deriv_3D","text":"kernel_deriv_3D(kernel::WendlandC8, u::Real, h_inv::Real)\n\nEvaluate the derivative of the WendlandC8 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_1D-Union{Tuple{T}, Tuple{Cubic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_1D","text":"kernel_value_1D(kernel::Cubic{T}, u::Real, h_inv::Real) where T\n\nEvaluate cubic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_1D-Union{Tuple{T}, Tuple{Quintic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_1D","text":"kernel_value_1D(kernel::Quintic{T}, u::Real, h_inv::Real) where T\n\nEvaluate quintic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_1D-Union{Tuple{T}, Tuple{WendlandC2{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_1D","text":"kernel_value_1D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC2 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_1D-Union{Tuple{T}, Tuple{WendlandC4{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_1D","text":"kernel_value_1D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC4 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_1D-Union{Tuple{T}, Tuple{WendlandC6{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_1D","text":"kernel_value_1D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC6 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_2D-Union{Tuple{T}, Tuple{Cubic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_2D","text":"kernel_value_2D(kernel::Cubic{T}, u::Real, h_inv::Real) where T\n\nEvaluate cubic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_2D-Union{Tuple{T}, Tuple{Quintic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_2D","text":"kernel_value_2D(kernel::Quintic{T}, u::Real, h_inv::Real) where T\n\nEvaluate quintic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_2D-Union{Tuple{T}, Tuple{WendlandC2{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_2D","text":"kernel_value_2D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC2 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_2D-Union{Tuple{T}, Tuple{WendlandC4{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_2D","text":"(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC4 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_2D-Union{Tuple{T}, Tuple{WendlandC6{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_2D","text":"(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC6 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_2D-Union{Tuple{T}, Tuple{WendlandC8{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_2D","text":"kernel_value_2D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC6 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_3D-Union{Tuple{T}, Tuple{Cubic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_3D","text":"kernel_value_3D(kernel::Cubic{T}, u::Real, h_inv::Real) where T\n\nEvaluate cubic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_3D-Union{Tuple{T}, Tuple{Quintic{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_3D","text":"kernel_value_3D(kernel::Quintic{T}, u::Real, h_inv::Real) where T\n\nEvaluate quintic spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_3D-Union{Tuple{T}, Tuple{WendlandC2{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_3D","text":"kernel_value_3D(kernel::WendlandC2{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC2 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_3D-Union{Tuple{T}, Tuple{WendlandC4{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_3D","text":"kernel_value_3D(kernel::WendlandC4{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC4 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_3D-Union{Tuple{T}, Tuple{WendlandC6{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_3D","text":"kernel_value_3D(kernel::WendlandC6{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC6 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.kernel_value_3D-Union{Tuple{T}, Tuple{WendlandC8{T}, Real, Real}} where T","page":"API reference","title":"SPHKernels.kernel_value_3D","text":"kernel_value_3D(kernel::WendlandC8{T}, u::Real, h_inv::Real) where T\n\nEvaluate WendlandC8 spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.δρ₁-Tuple{SPHKernel, Real, Real, Real}","page":"API reference","title":"SPHKernels.δρ₁","text":"δρ₁(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)\n\nCorrects the 1D density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.δρ₂-Tuple{SPHKernel, Real, Real, Real}","page":"API reference","title":"SPHKernels.δρ₂","text":"δρ₂(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)\n\nCorrects the 2D density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.δρ₃-Tuple{SPHKernel, Real, Real, Real}","page":"API reference","title":"SPHKernels.δρ₃","text":"δρ₃(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)\n\nCorrects the 3D density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.∇𝒲₁-Tuple{SPHKernel, Real, Real}","page":"API reference","title":"SPHKernels.∇𝒲₁","text":"∇𝒲₁( kernel::SPHKernel, u::Real, h_inv::Real)\n\nEvaluate 1D derivative at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.∇𝒲₂-Tuple{SPHKernel, Real, Real}","page":"API reference","title":"SPHKernels.∇𝒲₂","text":"∇𝒲₂( kernel::SPHKernel, u::Real, h_inv::Real)\n\nEvaluate 1D derivative at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.∇𝒲₃-Tuple{SPHKernel, Real, Real}","page":"API reference","title":"SPHKernels.∇𝒲₃","text":"∇𝒲₃( kernel::SPHKernel, u::Real, h_inv::Real)\n\nEvaluate 1D derivative at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.𝒲₁-Tuple{SPHKernel, Real, Real}","page":"API reference","title":"SPHKernels.𝒲₁","text":"𝒲₁( kernel::SPHKernel, u::Real, h_inv::Real)\n\nEvaluate 1D spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.𝒲₂-Tuple{SPHKernel, Real, Real}","page":"API reference","title":"SPHKernels.𝒲₂","text":"𝒲₂( kernel::SPHKernel, u::Real, h_inv::Real)\n\nEvaluate 2D spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#SPHKernels.𝒲₃-Tuple{SPHKernel, Real, Real}","page":"API reference","title":"SPHKernels.𝒲₃","text":"𝒲₃( kernel::SPHKernel, u::Real, h_inv::Real)\n\nEvaluate 3D spline at position u = fracxh.\n\n\n\n\n\n","category":"method"},{"location":"api/#Private-Functions","page":"API reference","title":"Private Functions","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [SPHKernels]\nPublic = false\nOrder = [:function]","category":"page"},{"location":"api/#Private-Types","page":"API reference","title":"Private Types","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [SPHKernels]\nPublic = false\nOrder = [:type]","category":"page"},{"location":"extending/#Adding-Kernels","page":"Extending","title":"Adding Kernels","text":"","category":"section"},{"location":"extending/","page":"Extending","title":"Extending","text":"CurrentModule = SPHKernels\nDocTestSetup = quote\n    using SPHKernels\nend","category":"page"},{"location":"extending/","page":"Extending","title":"Extending","text":"If you need a different kernel function than the ones I implemented you can add them by defining a new kernel struct as a subtype of SPHKernel","category":"page"},{"location":"extending/","page":"Extending","title":"Extending","text":"using SPHKernels # hide\nstruct MyKernel{T} <: SPHKernel\n    n_neighbours::Int64\n    norm_1D::T\n    norm_2D::T\n    norm_3D::T\nend\n\n\"\"\"\n    MyKernel(T::DataType=Float64, n_neighbours::Integer=295)\n\nSet up a `MyKernel` kernel for a given DataType `T`.\n\"\"\"\nMyKernel(T::DataType=Float64, n_neighbours::Integer=1000) = MyKernel{T}(n_neighbours, T(1), T(1), T(1))","category":"page"},{"location":"extending/","page":"Extending","title":"Extending","text":"and defining its value and derivative in 2D and 3D, e.g.","category":"page"},{"location":"extending/","page":"Extending","title":"Extending","text":"function kernel_value_2D(kernel::MyKernel{T}, u::Real, h_inv::Real) where T\n\n    if u < 1\n        n = kernel.norm_2D * h_inv^2\n        return 1.0 * n |> T\n    else\n        return 0.0 |> T\n    end\n\nend","category":"page"},{"location":"extending/","page":"Extending","title":"Extending","text":"@inline function kernel_deriv_2D(kernel::MyKernel{T}, u::Real, h_inv::Real) where T\n    return 0.0 |> T\nend","category":"page"},{"location":"extending/","page":"Extending","title":"Extending","text":"so you can run the following code","category":"page"},{"location":"extending/","page":"Extending","title":"Extending","text":"k = MyKernel()\nu = 0.5\nh_inv = 1.0\n\nv = kernel_value_2D(k, u, h_inv)\nprintln(\"MyKernel value: $v\")\n\nd = kernel_deriv_2D(k, u, h_inv)\nprintln(\"MyKernel derivative: $d\")","category":"page"},{"location":"extending/#Contributing","page":"Extending","title":"Contributing","text":"","category":"section"},{"location":"extending/","page":"Extending","title":"Extending","text":"Please feel free to create a pull request if you feel that your kernel could be useful to others! I only ask you to also add unit tests for your kernels in the test/runtests.jl file.","category":"page"},{"location":"install/#Install","page":"Install","title":"Install","text":"","category":"section"},{"location":"install/","page":"Install","title":"Install","text":"As usual with Julia you can install the package via the internal package manager","category":"page"},{"location":"install/","page":"Install","title":"Install","text":"julia> ]\npkg> add https://github.com/LudwigBoess/SPHKernels.jl","category":"page"},{"location":"#Table-of-contents","page":"Table of Contents","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Table of Contents","title":"Table of Contents","text":"CurrentModule = SPHKernels\nDocTestSetup = quote\n    using SPHKernels\nend","category":"page"},{"location":"","page":"Table of Contents","title":"Table of Contents","text":"Pages = [ \"index.md\", \n          \"install.md\", \n          \"kernels.md\",\n          \"extending.md\", \n          \"api.md\"]\nDepth = 3","category":"page"}]
}
