"""
    This package contains a number of kernels frequently used in SPH, as well as their functions and derivatives an 2D and 3D.
"""

module SPHKernels

    export  kernel_value_1D,     𝒲₁,
            kernel_value_2D,     𝒲₂,
            kernel_value_3D,     𝒲₃,
            kernel_deriv_1D,    d𝒲₁,
            kernel_deriv_2D,    d𝒲₂,
            kernel_deriv_3D,    d𝒲₃,
            bias_correction_1D, δρ₁,
            bias_correction_2D, δρ₂,
            bias_correction_3D, δρ₃,
            SPHKernel,
            Cubic, 
            Quintic,
            WendlandC2,
            WendlandC4,
            WendlandC6,
            WendlandC8
        

    """
        SPHKernel

    Supertype for all SPH kernels.
    """
    abstract type SPHKernel end

    include("bsplines/Cubic.jl")
    include("bsplines/Quintic.jl")
    include("wendland/C2.jl")
    include("wendland/C4.jl")
    include("wendland/C6.jl")
    include("wendland/C8.jl")

    # multiple dispatch for nicer look

    """
        𝒲₁( kernel::SPHKernel, u::Real, h_inv::Real)

    Evaluate 1D spline at position ``u = \\frac{x}{h}``.
    """
    𝒲₁( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value_1D(kernel, u, h_inv)
    
    """
        d𝒲₁( kernel::SPHKernel, u::Real, h_inv::Real)

    Evaluate 1D derivative at position ``u = \\frac{x}{h}``.
    """
    d𝒲₁(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv_1D(kernel, u, h_inv)

    """
        𝒲₂( kernel::SPHKernel, u::Real, h_inv::Real)

    Evaluate 2D spline at position ``u = \\frac{x}{h}``.
    """
    𝒲₂( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value_2D(kernel, u, h_inv)

    """
        d𝒲₂( kernel::SPHKernel, u::Real, h_inv::Real)

    Evaluate 1D derivative at position ``u = \\frac{x}{h}``.
    """
    d𝒲₂(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv_2D(kernel, u, h_inv)

    """
        𝒲₃( kernel::SPHKernel, u::Real, h_inv::Real)

    Evaluate 3D spline at position ``u = \\frac{x}{h}``.
    """
    𝒲₃( kernel::SPHKernel, u::Real, h_inv::Real) = kernel_value_3D(kernel, u, h_inv)
    
    """
        d𝒲₃( kernel::SPHKernel, u::Real, h_inv::Real)

    Evaluate 1D derivative at position ``u = \\frac{x}{h}``.
    """
    d𝒲₃(kernel::SPHKernel, u::Real, h_inv::Real) = kernel_deriv_3D(kernel, u, h_inv)

    """ 
        δρ₁(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)

    Corrects the 1D density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
    """
    δρ₁(kernel::SPHKernel, density::Real, m::Real, h_inv::Real) = bias_correction_1D(kernel, density, m, h_inv)

    """ 
        δρ₂(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)

    Corrects the 2D density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
    """
    δρ₂(kernel::SPHKernel, density::Real, m::Real, h_inv::Real) = bias_correction_2D(kernel, density, m, h_inv)

    """ 
        δρ₃(kernel::SPHKernel, density::Real, m::Real, h_inv::Real)

    Corrects the 3D density estimate for the kernel bias. See Dehnen&Aly 2012, eq. 18+19.
    """
    δρ₃(kernel::SPHKernel, density::Real, m::Real, h_inv::Real) = bias_correction_3D(kernel, density, m, h_inv)
    
end # module
