"""
    kernel_norm(kernel::WendlandC2_1D{T}, h_inv::Real) where {T}

Calculate the normalisation factor for the WendlandC2 kernel.
"""
kernel_norm(kernel::Union{WendlandC2_1D{T}, WendlandC4_1D{T}, WendlandC6_1D{T}, WendlandC8_1D{T}}, 
            h_inv::Real) where {T} = kernel.norm * h_inv

            
"""
    kernel_norm(kernel::WendlandKernel, h_inv::Real) where {T}

Calculate the normalisation factor for the WendlandC2 kernel.
"""
function kernel_norm(kernel::WendlandKernel, h_inv::Real)
    if kernel.dim == Int8(2)
        return kernel.norm * h_inv*h_inv
    else
        return kernel.norm * h_inv*h_inv*h_inv
    end
end