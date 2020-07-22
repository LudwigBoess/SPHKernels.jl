"""
    This package contains a number of kernels frequently used in SPH, as well as their functions and derivatives an 2D and 3D.
"""

module SPHKernels

    export  kernel_value_2D,
            kernel_value_3D,
            kernel_deriv_2D,
            kernel_deriv_3D,
            bias_correction_2D,
            bias_correction_3D,
            SPHKernel,
            Cubic, 
            Quintic,
            WendlandC4,
            WendlandC6
        

    """
        SPHKernel

    Supertype for all SPH kernels.
    """
    abstract type SPHKernel end

    include("bsplines.jl")
    include("wendland.jl")

end # module
